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

rule fc_summary:
    input:
        "logs/quant/featurecounts.log"
    output:
        "04.Quant_featureCounts/fc.summary"
    shell:
    """
echo "sample\tassignrate" > {output}
sed ':a;N;$!ba;s/\n//g' {input}| sed -e 's/Process BAM file /\n/g'|sed 1d|sed -e 's/.bam.*(/\t/g' -e 's/%.*//g' >> {output}
    """