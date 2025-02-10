rule fastp_qc:
    input:
        r1 = "01.Rawdata/{sample}/*_1.fq.gz",
        r2 = "01.Rawdata/{sample}/*_2.fq.gz"
    output:
        r1_clean = "02.QC_fastp/{sample}/{sample}_R1.clean.fastq.gz",
        r2_clean = "02.QC_fastp/{sample}/{sample}_R2.clean.fastq.gz",
        html = "02.QC_fastp/reports/{sample}.html",
        json = "02.QC_fastp/reports/{sample}.json",
    log:
        "logs/qc/{sample}.log"
    threads: 2
    resources:
        mem_mb=2048,
        runtime=60
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
        "02.QC_fastp/QC_summary_table.tsv",
        "02.QC_fastp/multiqc_report.html"
    shell:
        """
        # Generate tabular summary
        module load GCC/12.2.0  OpenMPI/4.1.4 MultiQC/1.14

        echo -e "Sample\tTotal_Reads\tQ20_R1\tQ20_R2\tGC_R1\tGC_R2" > {output[0]}
        for f in {input}; do
            sample=$(basename $f .json)
            total_reads=$(jq '.summary.before_filtering.total_reads' $f)
            q20_r1=$(jq '.summary.before_filtering.q20_rate' $f)
            q20_r2=$(jq '.summary.before_filtering.q20_rate' $f)
            gc_r1=$(jq '.summary.before_filtering.gc_content' $f)
            gc_r2=$(jq '.summary.before_filtering.gc_content' $f)
            echo -e "$sample\t$total_reads\t$q20_r1\t$q20_r2\t$gc_r1\t$gc_r2" >> {output[0]}
        done
        
        # Generate MultiQC report
        multiqc 02.QC/reports/ -o 02.QC/
        """

