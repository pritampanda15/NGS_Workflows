# Assessing heterozygosity of species, if DNA-seq exist
if DATA_DICT['DNA']:
    # k-mer histogram
    rule jellyfish_count_dna:
        input:
            expand(opj(OUTDIR, "trimmed/{dna}/{dna}.clean.{run}.fq.gz"), dna=DATA_DICT['DNA'], run=RUN),
        output:
            opj(OUTDIR, "jellyfish/DNA.jf"),
        log:
            opj(OUTDIR, "logs/jellyfish/count_DNA.log")
        params:
            kmer = 31,
            extra = "--canonical",  # You should always use "canonical k-mers" (-C) since the sequencing reads will come from both the forward and reverse strand of DNA.
        threads: 10
        resources:
            mem_gb = '1G',
        wrapper:
            get_wrapper('jellyfish', 'count')


# HIFI k-mer histogram
rule jellyfish_count_hifi:
    input:
        expand(opj(OUTDIR, "trimmed/{hifi}/{hifi}.filt.fastq.gz"), hifi=DATA_DICT['HIFI']),
    output:
        opj(OUTDIR, "jellyfish/HIFI.jf"),
    log:
        opj(OUTDIR, "logs/jellyfish/count_HIFI.log")
    params:
        kmer = 31,
        extra = "--canonical",  # You should always use "canonical k-mers" (-C) since the sequencing reads will come from both the forward and reverse strand of DNA.
    threads: 10
    resources:
        mem_gb = '1G',
    wrapper:
        get_wrapper('jellyfish', 'count')


rule jellyfish_hist:
    input:
        opj(OUTDIR, "jellyfish/{access}.jf"),
    output:
        opj(OUTDIR, "jellyfish/{access}_jf.hist"),
    log:
        opj(OUTDIR, "logs/jellyfish/hist_{access}.log")
    threads: 2
    wrapper:
        get_wrapper('jellyfish', 'histo')

# plot hist
rule genomescope:
    input:
        rules.jellyfish_hist.output
    output:
        directory(opj(OUTDIR, "genomescope/{access}/v1"))
    log:
        opj(OUTDIR, "logs/genomescope/{access}_v1.log")
    params:
        kmer = 31,
        read_length = 150,
        kmer_max = -1,
        verbose = True,
    wrapper:
        get_wrapper('scripts', 'R', 'genomescope')

rule genomescope2:
    input:
        rules.jellyfish_hist.output
    output:
        directory(opj(OUTDIR, "genomescope/{access}/v2"))
    log:
        opj(OUTDIR, "logs/genomescope/{access}_v2.log")
    params:
        kmer = 31,
        extra = ""
    wrapper:
        get_wrapper('genomescope2')