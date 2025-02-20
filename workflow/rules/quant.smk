localrules: fc_merge, fc_summary, fc_plot

rule featurecounts:
    input:
        "03.Alignment_hisat2/{sample}/{sample}.bam"
    output:
        rc = "04.Quant_featureCounts/{sample}/{sample}.counts_raw.tsv",
        summ = "04.Quant_featureCounts/{sample}/{sample}.counts_raw.tsv.summary"
    log:
        "logs/quant/{sample}.featurecounts.log"
    params:
        anno = config['genome']['annotation']
    resources:
        runtime = 60,
        cpus_per_task = 4,
        mem_mb = 4096
    shell:
        """
module load GCC/12.3.0 Subread/2.0.8

featureCounts \
-T {resources.cpus_per_task} \
-a {params.anno} \
-o {output.summ} \
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

rule fc_merge:
    input:
        expand("04.Quant_featureCounts/{sample}/{sample}.counts_raw.tsv", sample=config['samples']),
    output:
        "04.Quant_featureCounts/counts_raw.tsv",
        "04.Quant_featureCounts/counts_filter.tsv"
    params:
        pipepath = config['pipepath']
    shell:
        """
echo {input}
echo {output}
{params.pipepath}/scripts/raw.rc.filter.sh {input} {output}
        """

rule fc_summary:
    input:
        expand("logs/quant/{sample}.featurecounts.log", sample=config['samples'])
    output:
        "04.Quant_featureCounts/fc.summary"
    shell:
        """
echo {input}
echo "sample\tassignrate" > {output}
sed ':a;N;$!ba;s/\\n//g' {input}| sed -e 's/Process BAM file /\\n/g'|sed 1d|sed -e 's/.bam.*(/\\t/g' -e 's/%.*//g' >> {output}
        """

rule fc_plot:
    input:
        expand("04.Quant_featureCounts/{sample}/{sample}.counts_raw.tsv.summary", sample=config['samples'])
    output:
        "04.Quant_featureCounts/assignment_stacked_barplot.png"
    params:
        pipepath = config['pipepath']
    shell:
        """
module load GCC/12.2.0 OpenMPI/4.1.4 R/4.3.1
export R_LIBS_USER="/scratch/group/lilab/software/R_library/4.3"
Rscript {params.pipepath}/scripts/fc.assign.plot.R {input} {output}
        """

