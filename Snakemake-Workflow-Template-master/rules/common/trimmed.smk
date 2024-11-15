if TRIMMING == "fastp":
    rule fastp:
        input:
            reads = get_fastq,
        output:
            trimmed = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
            failed = opj(OUTDIR, "trimmed/{sample}/{sample}.failed.fastq"),
            html = opj(OUTDIR, "trimmed/report/{sample}.fastp.html"),
            json = opj(OUTDIR, "trimmed/report/{sample}.fastp.json")
        log:
            opj(OUTDIR, "logs/trimmed/fastp_{sample}.log")
        params:
            extra = config['parameters']['fastp']['extra'],  # optional parameters
            adapters = get_adapter('fastp')
        threads: 2
        wrapper:
            get_wrapper('fastp')
elif TRIMMING == "trimmomatic":
    rule trimmomatic:
        input:
            get_fastq,
        output:
            trimmed = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
            # reads where trimming entirely removed the mate
            unpaired = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.unpaired.{run}.fq.gz"), run=RUN),
        log:
            opj(OUTDIR, "logs/trimmed/trimmomatic_{sample}.log")
        params:
            # list of trimmers (see manual)
            trimmer = config['parameters']['trimmomatic']['trimmer'],
            # get adapters file
            extra = get_adapter('trimmomatic'),
            compression_level = config['parameters']['trimmomatic']['compression_level']
        threads:
            32
        resources:
            mem_mb = config['parameters']['trimmomatic']['mem']
        wrapper:
            get_wrapper("trimmomatic")
elif TRIMMING == "cutadapt":
    rule cutadapt:
        input:
            get_fastq,
        output:
            trimmed = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
            qc = opj(OUTDIR, "trimmed/report/{sample}.qc.txt"),
        params:
            # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
            adapters = get_adapter('cutadapt'),
            # https://cutadapt.readthedocs.io/en/stable/guide.html#
            extra = config['parameters']['cutadapt']['extra'],
        log:
            opj(OUTDIR, "logs/trimmed/cutadapt_{sample}.log"),
        threads: 4  # set desired number of threads here
        wrapper:
            get_wrapper("cutadapt")
elif TRIMMING == "trim_galore":
    rule trim_galore:
        input:
            get_fastq,
        output:
            trimmed = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
            report = expand(opj(OUTDIR, "trimmed/report/{{sample}}.{run}.trimming_report.txt"), run=RUN),
        params:
            extra = config['parameters']['trim_galore']['extra'],
        log:
            opj(OUTDIR, "logs/trimmed/trim_galore_{sample}.log"),
        wrapper:
            get_wrapper("trim_galore")
else:
    raise ValueError(f'the rule: {TRIMMING} not support!')