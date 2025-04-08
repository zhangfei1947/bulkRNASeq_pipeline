localrules: fc_summary, fc_plot

rule featurecounts:
    input:
        expand("03.Alignment_hisat2/{sample}/{sample}.bam", sample=config['samples'])
    output:
        "04.Quant_featureCounts/counts_raw.tsv",
        "04.Quant_featureCounts/counts_raw.tsv.summary"
    log:
        "logs/quant/featurecounts.log"
    params:
        anno = config['genome']['annotation']
    resources:
        runtime = 360,
        cpus_per_task = int(len(config['samples'])/3),
        mem_mb = 400*len(config['samples'])
    shell:
        """
featureCounts \
-T {resources.cpus_per_task} \
-a {params.anno} \
-o {output[0]} \
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

rule fc_filter:
    input:
        "04.Quant_featureCounts/counts_raw.tsv"
    output:
        "04.Quant_featureCounts/counts_filter.tsv"
    params:
        pipepath = config['pipepath']
    container: None
    script:
        "scripts/raw.rc.filter.sh {input} {output}"

rule fc_summary:
    input:
        "logs/quant/featurecounts.log"
    output:
        "04.Quant_featureCounts/fc.summary"
    container: None
    shell:
        """
echo "sample\tassignrate" > {output}
sed ':a;N;$!ba;s/\\n//g' {input}| sed -e 's/Process BAM file /\\n/g'|sed 1d|sed -e 's/.bam.*(/\\t/g' -e 's/%.*//g' >> {output}
        """

rule fc_plot:
    input:
        "04.Quant_featureCounts/counts_raw.tsv.summary"
    output:
        "04.Quant_featureCounts/assignment_stacked_barplot.png"
    params:
        pipepath = config['pipepath']
    script:
        "scripts/fc.assign.plot.R {input} {output}"

