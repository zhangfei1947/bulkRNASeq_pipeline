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
        runtime = 30,
        nodes = 1,
        ntasks = 1,
        cpus_per_task= 2,
        mem_mb = 4000
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
        # Generate tabular summary
        alias jq="/scratch/group/lilab/software/jq-linux-i386"
        echo {input}
        ../scripts/jq.qc.sum.sh {output}
        """

rule duprate_plot:
    input:
        "02.QC_fastp/QC_summary_table.tsv"
    output:
        "02.QC_fastp/duprate.boxplot.png"
    shell:
        """
        module load GCC/12.2.0 OpenMPI/4.1.4 R/4.3.1
        export R_LIBS_USER="/scratch/group/lilab/software/R_library/4.3"
        Rscript ../scripts/duprate.plot.R {input}
        """




