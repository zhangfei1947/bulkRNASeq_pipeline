
rule featurecounts:
    input:
        bams = expand("03.Alignment_hisat2/{sample}/{sample}.bam", sample=config['samples']),
        anno = config['genome']['annotation']
    output:
        "04.Quant_featureCounts/counts_raw.tsv"
    log:
        "logs/quant/featurecounts.log"
    threads: 2
    params:
        strand = 2
    shell:
        """
        module load GCC/12.3.0 Subread/2.0.8

        featureCounts \
        -T {threads} \
        -a {input.anno} \
        -o {output} \
        -F GTF -t exon -g gene_id \
        -s {params.strand} \
        -p \
        --countReadPairs \
        -O -M --fraction \
        -P -B \
        -d 40 \
        -D 800 \
        {input.bams} > {log} 2>&1
        """
