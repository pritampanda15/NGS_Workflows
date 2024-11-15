rule bamCoverage:
    input:
        # Required input, must index
        rules.downsample.output if DOWNSAMPLE else rules.filterChrM.output,
    output:
        # Required output.
        # Output file format should be one of ['bw', 'bigwig', 'bigWig', 'bedgraph', 'bedGraph'].
        opj(OUTDIR, "macs2/bigwig/{sample}.norm.bw")
    params:
        # Optional parameters.
        extra = config['parameters']['deeptools']['bamcoverage']['extra'] + f' --effectiveGenomeSize {CHROM_SIZE}',
        # extar = '--binSize 1 --normalizeUsing RPKM',
    threads: 1
    log: 
        opj(OUTDIR, "logs/deeptools/bamcoverage_{sample}.log")
    wrapper:
        get_wrapper('deeptools', 'bamcoverage')

rule TSSEnrichment:
    input: 
        # bigwig = "macs2/bigwig/{pair}.norm.bw", 
        bigwig = get_bigwig,
        bed = config['data']['gtf']
    output:
        matrix_gz = opj(OUTDIR, "macs2/matrix/{control}.tss.matrix.gz"),   # required
    log:
        opj(OUTDIR, "logs/deeptools/compute_tss_matrix_{control}.log")
    params:
        # required argument, choose "scale-regions" or "reference-point"
        command = "reference-point",
        # optional parameters
        extra = config['parameters']['deeptools']['computematrix']['tss'],
        labels = lambda wildcards, input: [f.split('/')[-1].strip('.norm.bw') for f in input.bigwig]
    wrapper:
        get_wrapper("deeptools", "computematrix")

rule plotTSSHeatmap:
    input:
         rules.TSSEnrichment.output.matrix_gz
    output:
        heatmap_img = opj(OUTDIR, "macs2/matrix/{control}.tss.heatmap.pdf"),  # required
    log:
        opj(OUTDIR, "logs/deeptools/TSS_heatmap_{control}.log")
    params:
        # optional parameters
        extra = config['parameters']['deeptools']['plotheatmap'],
    wrapper:
        get_wrapper("deeptools", "plotheatmap")

rule plotTSSProfile:
    input:
        rules.TSSEnrichment.output.matrix_gz
    output:
        plot_img = opj(OUTDIR, "macs2/matrix/{control}.tss.profile.pdf"),  # required
    log:
        opj(OUTDIR, "logs/deeptools/TSS_profile_{control}.log")
    params:
        # optional parameters
        extra = config['parameters']['deeptools']['plotprofile'],
    wrapper:
        get_wrapper("deeptools", "plotprofile")

rule genbodyEnrichment:
    input:
        bigwig = get_bigwig,
        bed = config['data']['gtf']
    output:
        matrix_gz = opj(OUTDIR, "macs2/matrix/{control}.genebody.matrix.gz"),   # required
    log:
        opj(OUTDIR, "logs/deeptools/compute_genebody_matrix_{control}.log")
    params:
        # required argument, choose "scale-regions" or "reference-point"
        command = "scale-regions",
        # optional parameters
        extra = config['parameters']['deeptools']['computematrix']['genbody'],
        labels = lambda wildcards, input: [f.split('/')[-1].strip('.norm.bw') for f in input.bigwig]
    wrapper:
        get_wrapper("deeptools", "computematrix")

rule plotBodyHeatmap:
    input:
         rules.genbodyEnrichment.output.matrix_gz
    output:
        heatmap_img = opj(OUTDIR, "macs2/matrix/{control}.genebody.heatmap.pdf"),  # required
    log:
        opj(OUTDIR, "logs/deeptools/genebody_heatmap_{control}.log")
    params:
        # optional parameters
        extra = config['parameters']['deeptools']['plotheatmap'],
    wrapper:
        get_wrapper("deeptools", "plotheatmap")

rule ploGenebodytProfile:
    input:
        rules.genbodyEnrichment.output.matrix_gz
    output:
        plot_img = opj(OUTDIR, "macs2/matrix/{control}.genebody.profile.pdf"),  # required
    log:
        opj(OUTDIR, "logs/deeptools/genebody_profile_{control}.log")
    params:
        # optional parameters
        extra = config['parameters']['deeptools']['plotprofile'],
    wrapper:
        get_wrapper("deeptools", "plotprofile")


if len(SAMPLES) > 1:
    rule multibigwigsummary:
        input:
            # Required input.
            # bigwig = get_bigwig,
            bigwig = expand(opj(OUTDIR, "macs2/bigwig/{sample}.norm.bw"), sample=PAIRS),
        output:
            # Required output.
            out = opj(OUTDIR, "macs2/bigwig/scores_per_bin.npz"),
            # raw_count = ''
        params:
            # Optional parameters.
            subcommand = 'bins',
            extra = config['parameters']['deeptools']['multibigwigsummary'].get('extra'),
            labels = lambda wildcards, input: [f.split('/')[-1].strip('.norm.bw') for f in input.bigwig]
        threads: 1
        log: 
            opj(OUTDIR, "logs/deeptools/multibigwigsummary.log")
        wrapper: 
            get_wrapper("deeptools", "multibigwigsummary")
            
    rule plotcorrelation:
        input:
            # Required input.
            rules.multibigwigsummary.output.out,
        output:
            img = opj(OUTDIR, "macs2/bigwig/heatmap_spearman_corr_readCounts.pdf"),
            # matrix = 'SpearmanCorr_readCounts.tab'
        params:
            # Optional parameters.
            extra = config['parameters']['deeptools']['plotcorrelation'].get('extra'),
            cor = config['parameters']['deeptools']['plotcorrelation'].get('cor'),
            title = config['parameters']['deeptools']['plotcorrelation'].get('title'),
        threads: 1
        log: 
            opj(OUTDIR, "logs/deeptools/plotcorrelation.log")
        wrapper: 
            get_wrapper("deeptools", "plotcorrelation")