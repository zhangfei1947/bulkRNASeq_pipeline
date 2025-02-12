localrules: fc_summary

rule featurecounts:
    input:
        expand("03.Alignment_hisat2/{sample}/{sample}.bam", sample=config['samples'])
    output:
        "04.Quant_featureCounts/counts_raw.tsv"
    log:
        "logs/quant/featurecounts.log"
    params:
        anno = config['genome']['annotation']
    threads: 4
    resources:
        runtime = 60,
        nodes = 1,
        ntasks = 1,
        cpus_per_task= 4,
        mem_mb = 6000
    shell:
        """
        module load GCC/12.3.0 Subread/2.0.8

        featureCounts \
        -T {threads} \
        -a {params.anno} \
        -o {output} \
        -F GTF -t exon -g gene_id \
        -s 2 \
        -p \
        --countReadPairs \
        -O -M --fraction \
        -P -B \
        -d 40 \
        -D 800 \
        {input} > {log} 2>&1
        """
