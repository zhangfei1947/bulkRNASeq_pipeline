localrules: qc_summary, duprate_plot

rule fastp_qc:
    input:
        r1 = "01.Rawdata/{sample}/{sample}_R1.fq.gz",
        r2 = "01.Rawdata/{sample}/{sample}_R2.fq.gz"
    output:
        r1_clean = "02.QC_fastp/{sample}/{sample}_R1.clean.fastq.gz",
        r2_clean = "02.QC_fastp/{sample}/{sample}_R2.clean.fastq.gz",
        html = "02.QC_fastp/reports/{sample}.html",
        json = "02.QC_fastp/reports/{sample}.json",
    log:
        "logs/qc/{sample}.log"
    threads: 2
    resources:
        cpus_per_task= 2,
        mem_mb = 4096
    shell:
        """
        module load GCC/11.2.0 fastp/0.23.2
        fastp --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1_clean} --out2 {output.r2_clean} \
            -j {output.json} -h {output.html} \
            --thread {threads} > {log} 2>&1
        """

rule qc_summary:
    input:
        expand("02.QC_fastp/reports/{sample}.json", sample=config['samples'])
    output:
        "02.QC_fastp/QC_summary_table.tsv"
    shell:
        """
export PATH="/scratch/group/lilab/software/jq:$PATH"

echo -e "Sample\tReads_brfore\tReads_after\tQ20_R1_before\tQ20_R2_before\tQ20_R1_after\tQ20_R2_after\tduplication_rate\tinsert_size" > {output}
for f in 02.QC_fastp/reports/*json; do
    sample=$(basename $f .json)
    total_reads_b=$(jq '.summary.before_filtering.total_reads' $f)
    total_reads_a=$(jq '.summary.after_filtering.total_reads' $f)
    q20b_r1=$(jq '.summary.before_filtering.q20_rate' $f)
    q20b_r2=$(jq '.summary.before_filtering.q20_rate' $f)
    q20a_r1=$(jq '.summary.after_filtering.q20_rate' $f)
    q20a_r2=$(jq '.summary.after_filtering.q20_rate' $f)
    dup_rate=$(jq '.duplication.rate' $f)
    insert_size=$(jq '.insert_size.peak' $f)
    echo -e "$sample\t$total_reads_b\t$total_reads_a\t$q20b_r1\t$q20b_r2\t$q20a_r1\t$q20a_r2\t$dup_rate\t$insert_size" >> {output}
done
        """

rule duprate_plot:
    input:
        "02.QC_fastp/QC_summary_table.tsv"
    output:
        "02.QC_fastp/duprate.boxplot.png"
    params:
        pipepath = config['pipepath']
    shell:
        """
module load GCC/12.2.0 OpenMPI/4.1.4 R/4.3.1
export R_LIBS_USER="/scratch/group/lilab/software/R_library/4.3"
Rscript {params.pipepath}/scripts/duprate.plot.R {input} {output}
        """




