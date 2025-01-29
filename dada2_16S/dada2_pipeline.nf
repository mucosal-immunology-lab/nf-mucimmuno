// dada2_pipeline.nf - Main Nextflow pipeline file

// Import the DEMULTIPLEX module
include { DEMULTIPLEX; USE_EXISTING_DEMULTIPLEXED } from './modules/demultiplex/demultiplex.nf'
include { QUALITY_CHECK; QUALITY_CHECK_AGGREGATE } from './modules/quality_check/quality_check.nf'
include { FILTER_AND_TRIM; FILTER_AND_TRIM_AGGREGATE } from './modules/filter_and_trim/filter_and_trim.nf'
include { ERROR_MODEL; ERROR_MODEL_AGGREGATE } from './modules/error_model/error_model.nf'
include { MERGE_PAIRED_ENDS; MERGE_SEQTABS; REMOVE_CHIMERAS; READ_TRACKING } from './modules/count_table_generation/count_table.nf'
include { ASSIGN_TAXONOMY; ASSIGN_SPECIES } from './modules/assign_taxonomy/assign_taxonomy.nf'
include { DE_NOVO_PHYLO_TREE } from './modules/phylo_tree/phylo_tree.nf'

// Define the input channel using the provided sample sheet
Channel.fromPath(params.sample_sheet)
    .splitCsv(header: true)
    .map { row -> tuple(row.run_name, file(row.folder_path)) }  // Ensures folder paths are correctly interpreted
    .set { input_samples }

// Print input_samples to the console for debugging
// input_samples.view()

// Define variables for each run based on folder contents
input_samples.map { run_name, folder_path ->
    def demultiplexed_folder = folder_path.resolve("demultiplexed")

    if (demultiplexed_folder.exists()) {
        // Pre-existing demultiplexed folder
        println "Found 'demultiplexed' folder for ${run_name}. Skipping DEMULTIPLEX."
        [run_name, demultiplexed_folder, true]  // Tag as pre-demultiplexed
    } else {
        // No demultiplexed folder; prepare for DEMULTIPLEX
        println "No 'demultiplexed' folder found for ${run_name}. Preparing for DEMULTIPLEX."
        def r1_file = folder_path.resolve("R1.fastq.gz")
        def r2_file = folder_path.resolve("R2.fastq.gz")
        def index_file = folder_path.resolve("Index.fastq.gz")
        def barcode_file = folder_path.listFiles()?.find { it.name.contains("barcode_to_sample") }
        def barcode_path = barcode_file ? file(barcode_file.toAbsolutePath()) : error("Barcode file not found in ${folder_path}.")
        [run_name, r1_file, r2_file, index_file, barcode_path, false]  // Tag as requiring DEMULTIPLEX
    }
}.set { run_files }

// Print run_files to the console for debugging
// run_files.view()

workflow {
    // ===============================
    // D E M U L T I P L E X   R U N S
    // ===============================

    // Split runs into two channels
    needs_demux = run_files.filter { it[5] == false }
    prelinked_runs = run_files.filter { it[2] == true }

    // Run DEMULTIPLEX for raw runs
    demultiplexed_results = needs_demux.map { run_name, r1_file, r2_file, index_file, barcode_file, _ ->
        tuple(run_name, r1_file, r2_file, index_file, barcode_file)
    } | DEMULTIPLEX

    // Link pre-existing demultiplexed runs
    prelinked_results = prelinked_runs.map { run_name, demultiplexed_folder, _ ->
        println "Linking pre-existing demultiplexed data for ${run_name}."
        tuple(run_name, demultiplexed_folder)
    }

    // Combine outputs into a unified channel
    all_demultiplexed = demultiplexed_results.mix(prelinked_results)

    // =========================================
    // P L O T   Q U A L I T Y   P R O F I L E S
    // =========================================

    // Use the demultiplexed output as input for QUALITY_CHECK
    quality_rds_files_channel = all_demultiplexed | QUALITY_CHECK

    // Collect all RDS files into a list
    quality_rds_files_list = quality_rds_files_channel.collect()

    // Aggregate all quality RDS files into a single PDF
    QUALITY_CHECK_AGGREGATE(quality_rds_files_list)

    // =================================================
    // F I L T E R   A N D   T R I M   S E Q U E N C E S
    // =================================================

    // Filter and trim the sequences
    filtered_files_channel = all_demultiplexed | FILTER_AND_TRIM
    filter_plot_rds_files_list = filtered_files_channel.filter_plots.collect()

    // Aggregate all filtering RDS files into a single PDF
    FILTER_AND_TRIM_AGGREGATE(filter_plot_rds_files_list)

    // ===============================================================
    // G E N E R A T E   S E Q U E N C I N G   E R R O R   M O D E L S
    // ===============================================================

    // Generate sequencing error models
    filtered_runs = filtered_files_channel.filtered_files
    error_models_channel = filtered_runs | ERROR_MODEL
    error_plot_rds_files_list = error_models_channel.error_plots.collect()

    // Aggregate all error model RDS files into a single PDF
    ERROR_MODEL_AGGREGATE(error_plot_rds_files_list)

    // Add error model to the filtered runs channel
    error_models_list = error_models_channel.error_model

    // =========================================================
    // G E N E R A T E   S E Q T A B S   F O R   E A C H   R U N
    // =========================================================

    // Create a tuple input for seqtab generation
    // Tuple will contain: [ run_name, filtered_runs, error_runs ]
    seqtab_input = filtered_runs.join(error_models_list, by: 0)

    // Paired-ends (PE) merging
    seqtab_channel = MERGE_PAIRED_ENDS(seqtab_input)
    seqtab_individual = seqtab_channel
    .map { run_name, seqtab -> seqtab }
    .collect()
    
    // =============================================
    // M E R G E   A L L   R U N S   T O G E T H E R
    // =============================================

    // Merge the individual seqtab data together
    merged_seqtab_channel = MERGE_SEQTABS(seqtab_individual)

    // =================================================
    // R E M O V E   C H I M E R I C   S E Q U E N C E S
    // =================================================
    
    // Remove chimeric sequences
    seqtab_no_chimeras_channel = REMOVE_CHIMERAS(merged_seqtab_channel)

    // =========================================
    // T R A C K   R E A D S   R E T E N T I O N
    // =========================================

    // Prepare runs_info
    runs_info_channel = filtered_runs.map { run_name, filtered_dir ->
        def dir_str = filtered_dir.toString()
        [ run_name, dir_str, "${dir_str}/filtering_report.rds" ]
    }

    // Collect runs_info as a list of runs
    runs_info_list_channel = runs_info_channel
    .collect()
    .map { runs_info ->
        // Split the flattened list into sublists of size 3 (run_name, filtered_dir, filtering_report)
        runs_info.collate(3)
    }

    // Generate the reads tracking plot
    READ_TRACKING(merged_seqtab_channel, seqtab_no_chimeras_channel.seqtab_nochim, runs_info_list_channel.collect())

    // ===================================================
    // A S S I G N   T A X O N O M Y   W I T H   S I L V A
    // ===================================================

    // Assign taxonomy
    taxonomy_channel = ASSIGN_TAXONOMY(seqtab_no_chimeras_channel.seqtab_nochim)
    taxonomy_species_channel = ASSIGN_SPECIES(taxonomy_channel.taxonomy_rds)

    // =================================================================
    // P R E P A R E   D E   N O V O   P H Y L O G E N E T I C   T R E E
    // =================================================================

    // Retrieve necessary elements
    tree_channel = DE_NOVO_PHYLO_TREE(seqtab_no_chimeras_channel.seqtab_nochim, taxonomy_species_channel.taxonomy_species_rds)
}


