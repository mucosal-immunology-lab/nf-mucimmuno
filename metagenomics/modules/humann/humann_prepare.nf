process PREPARE_HUMANN_DATABASES {
    tag "prepare_humann_dbs"
    label 'process_high'

    publishDir "$params.outdir/humann_db", mode: 'copy'

    output:
      path "humann_db", emit: humann_db

    script:
    """
    set -euo pipefail

    # helper: download with speed-timeout + fallback
    download_with_fallback(){
      url="\$1"; out="\$2"
      prefix="http://huttenhower.sph.harvard.edu/humann_data"
      mirror_prefix="https://g-227ca.190ebd.75bc.data.globus.org/humann_data"
      mirror="\${url/\$prefix/\$mirror_prefix}"

      echo "Downloading \$url → \$out (timeout=60s idle)"
      curl -L --fail --speed-limit 1 --speed-time 60 --connect-timeout 15 \
           --retry 0 --retry-delay 10 \
           -o "\$out" "\$url" \
      || {
        echo "Primary stalled or failed, fallback to mirror: \$mirror"
        curl -L --fail --retry 3 --retry-delay 10 --connect-timeout 15 \
             -o "\$out" "\$mirror"
      }
    }

    # 1) prepare base folders
    mkdir -p humann_db/chocophlan humann_db/uniref humann_db/utility_mapping

    # 2) dynamically grab the URLs
    URL_CHOCO=\$(humann_databases --available \
      | grep '^chocophlan : full ' \
      | awk -F' = ' '{print \$2}')
    URL_UNI=\$(humann_databases --available \
      | grep '^uniref : uniref90_diamond ' \
      | awk -F' = ' '{print \$2}')
    URL_UTIL=\$(humann_databases --available \
      | grep '^utility_mapping : full ' \
      | awk -F' = ' '{print \$2}')

    # 3) download + extract each with fallback
    echo
    download_with_fallback "\$URL_CHOCO" "humann_db/chocophlan/\$(basename \$URL_CHOCO)"
    (
      cd humann_db/chocophlan
      tar -xzf \$(basename "\$URL_CHOCO")
      rm -f \$(basename "\$URL_CHOCO")
    )

    echo
    download_with_fallback "\$URL_UNI" "humann_db/uniref/\$(basename \$URL_UNI)"
    (
      cd humann_db/uniref
      tar -xzf \$(basename "\$URL_UNI")
      rm -f \$(basename "\$URL_UNI")
    )

    echo
    download_with_fallback "\$URL_UTIL" "humann_db/utility_mapping/\$(basename \$URL_UTIL)"
    (
      cd humann_db/utility_mapping
      tar -xzf \$(basename "\$URL_UTIL")
      rm -f \$(basename "\$URL_UTIL")
    )

    # 4) update HUMAnN3 config
    echo
    humann_config --update database_folders nucleotide      "\${PWD}/humann_db/chocophlan"
    humann_config --update database_folders protein         "\${PWD}/humann_db/uniref"
    humann_config --update database_folders utility_mapping "\${PWD}/humann_db/utility_mapping"

    # 5) final sanity check
    echo
    echo "Contents of humann_db/:"
    ls -R humann_db
    echo
    echo "HUMAnN3 will use:"
    humann_config --print | grep database_folders
    """
}

process UPDATE_HUMANN_CONFIG {
    tag "update_humann_config"
    label 'process_high'
    errorStrategy 'terminate'

    input:
        /// the base HUMAnN3 database folder provided by the user
        val(humann_db)

    output:
        /// re-emit the same path so downstream steps keep using it
        val(humann_db), emit: humann_db

    script:
    """
    set -euo pipefail

    echo "Updating HUMAnN3 config to use provided databases at: ${humann_db}"

    # nucleotide (ChocoPhlAn)
    humann_config --update database_folders nucleotide    "${humann_db}/chocophlan"

    # protein (UniRef90 DIAMOND)
    humann_config --update database_folders protein       "${humann_db}/uniref"

    # utility mapping
    humann_config --update database_folders utility_mapping "${humann_db}/utility_mapping"

    echo "HUMAnN3 configuration is now:"
    humann_config --print | grep database_folders
    """
}

process PREPARE_METAPHLAN_DB {
  tag "prepare_metaphlan_db"
  publishDir "$params.outdir/metaphlan_db", mode: 'copy'

  output:
  path "metaphlan_db", emit: metaphlan_db

  script:
  """
  set -euo pipefail
  mkdir -p metaphlan_db
  metaphlan --install --bowtie2db metaphlan_db
  """
}

process MERGE_READS {
    tag "${meta.id}"
    label 'process_high'
    
    input:
      tuple val(meta), path(reads)

    output:
      tuple val(meta), path("${meta.id}.merged.fastq.gz"), emit: merged_reads

    script:
    """
    set -euo pipefail

    echo "[${meta.id}] Merging reads: ${reads}"
    cat ${reads} > ${meta.id}.merged.fastq.gz

    # quick sanity check
    [ -s ${meta.id}.merged.fastq.gz ] || {
      echo "Error: merged FASTQ is empty!" >&2
      exit 1
    }
    """
}