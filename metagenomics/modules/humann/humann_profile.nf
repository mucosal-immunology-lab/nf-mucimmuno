process HUMANN_PROFILE {
    tag { meta.id }
    label 'process_high'
    conda "${moduleDir}/environment.yaml"

    input:
      tuple val(meta), path(merged_fastq)
      val(metaphlan_db)
      path(humann_db)

    output:
      // gene families
      tuple val(meta), path("${meta.id}_genefamilies.tsv"),    emit: genefam
      // pathway coverage
      tuple val(meta), path("${meta.id}_pathcoverage.tsv"),    emit: pathcov
      // pathway abundance
      tuple val(meta), path("${meta.id}_pathabundance.tsv"),   emit: pathabun

    script:
    """
    set -euo pipefail

    # Detect the one subdirectory under metaphlan_db that starts with "mpa_"
    MPA_IDX=\$(ls ${metaphlan_db} | grep '^mpa_' | head -n1)

    echo "→ Using Metaphlan index dir: \${MPA_IDX}"

    humann \\
      --input  ${merged_fastq} \\
      --output . \\
      --threads ${task.cpus} \\
      --output-basename ${meta.id} \\
      --bowtie-options="--very-sensitive-local" \\
      --metaphlan-options="--index ${metaphlan_db}/\${MPA_IDX}"
    """
}
