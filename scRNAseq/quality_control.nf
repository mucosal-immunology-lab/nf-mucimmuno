// IMPORT PROCESSES

include { SOUPX_CORRECTION; DECOMPRESS_EXPORT_ARCHIVE } from './modules/soupX_correction/main.nf'

workflow {
    // Define your input parameters
    export_archive = params.soupx_correction.export_archive ?: './results/export_files.tar.gz'
    data_option = params.soupx_correction.data_option ?: 'GeneFull_Ex50pAS'
    pca_dims = params.soupx_correction.pca_dims ?: '1:30'

    // ===========================================================================================
    // R U N    S O U P X    C O R R E C T I O N    T O    R E M O V E    A M B I E N T    M R N A
    // ===========================================================================================

    // Channel to hold the path to the compressed file
    export_file_channel = Channel.fromPath(export_archive)

    // Decompress the export directory
    decompressed_dir_channel = DECOMPRESS_EXPORT_ARCHIVE(export_file_channel)

    // Find all sample directories within the decompressed directory
    sample_dirs_channel = decompressed_dir_channel
        .map { file('results/export_directory').listFiles() }
        .flatten()
        .filter { it.isDirectory() }
        .view { "Sample directory: ${it}" }

    // Prepare parameters for each sample directory
    soupx_correction_channel = sample_dirs_channel
        .map { sample_dir -> 
            def sample_name = sample_dir.getName()
            tuple(sample_dir, data_option, pca_dims, [id: sample_name])
        }
        .view { "Sample params: ${it}" }

    // Execute the SOUPX_CORRECTION process
    soupX_folders = Channel.empty()
    ch_soupX = SOUPX_CORRECTION(soupx_correction_channel)
    soupX_folders = ch_soupX.soupX_folders
        .view()
}