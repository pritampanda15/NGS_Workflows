nextflow.enable.dsl=2

params.sra_ids = ["SRR7179504", "SRR7179505", "SRR7179506", "SRR7179507",
    "SRR7179508", "SRR7179509", "SRR7179510", "SRR7179511",
    "SRR7179520", "SRR7179521", "SRR7179522", "SRR7179523",
    "SRR7179524", "SRR7179525", "SRR7179526", "SRR7179527",
    "SRR7179536", "SRR7179537", "SRR7179540","SRR7179541"]
    
params.prefetch_results = "prefetch_results"
params.fastqdump_results = "fastqdump_results"

sra_ids_ch = Channel.fromList(params.sra_ids)

workflow {
    sra_ids_mapped = sra_ids_ch.map { it.toString() }
    sra_files = prefetchSRA(sra_ids_mapped)
    fastq_files = fastqDump(sra_files)
}

process prefetchSRA {
    publishDir "${params.prefetch_results}", mode: 'copy'
    cleanup = true
    input:
    val sra_id 

    output:
    path "results/sra/${sra_id}/${sra_id}.sra"

    script:
    """
    mkdir -p results/sra
    echo "Prefetching SRA ID: ${sra_id}"
    prefetch -O results/sra/ ${sra_id}
    if [ -f results/sra/${sra_id}/${sra_id}.sra ]; then
        echo "File downloaded successfully."
    else
        echo "Download failed or incomplete."
        exit 1
    fi
    """
}

process fastqDump {
publishDir "${params.fastqdump_results}", mode: 'copy'
    cleanup = true
    input:
    path sra_file 

    output:
    path "fastq/${sra_file.getBaseName()}_*.fastq.gz"

    script:
    """
    mkdir -p fastq
    sra_id=\$(basename ${sra_file} .sra)
    fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip ${sra_file}
    mv fastq/*.fastq.gz fastq/\${sra_id}_*.fastq.gz
    """
}


process organizeResults {
    input:
    path sra_files
    path fastq_files

    script:
    """
    mkdir -p ${params.prefetch_results}
    mkdir -p ${params.fastqdump_results}

    mv ${sra_files} ${params.prefetch_results}/
    mv ${fastq_files} ${params.fastqdump_results}/
    """
}

