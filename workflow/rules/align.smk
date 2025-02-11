localrules: hisat2_summary

rule hisat2_align:
    input:
        r1 = "02.QC_fastp/{sample}_R1.clean.fastq.gz",
        r2 = "02.QC_fastp/{sample}_R2.clean.fastq.gz"
    output:
        "03.Alignment_hisat2/{sample}/{sample}.bam",
        "03.Alignment_hisat2/{sample}/{sample}.summary",
        "03.Alignment_hisat2/{sample}/{sample}.bam.bai"
    log:
        "logs/align/{sample}.log"
    threads: 2
    resources:
        mem_mb = 4000,
        runtime = 60,
        nodes = 1,
        ntasks_per_node = 2
    params:
        index = config['genome']['index'],
        splicesites = config['genome']['splicesites'],
        extra = config['align']
    shell:
        """
        module load Anaconda3/2024.02-1
        source activate hisat2

        hisat2 -x {params.index} --known-splicesite-infile {params.splicesites} \
            -p {threads} \
            --rna-strandness FR \
            --summary-file {output[1]} \
            --no-unal \
            -1 {input.r1} -2 {input.r2} 2> {log} |
        samtools view -@ {threads} -Sb |
        samtools sort -@ {threads} -o {output[0]}
        samtools index {output[0]}
        """
#rule hisat2_summary:
#    input:
#        "03.Alignment_hisat2/{sample}/{sample}.summary"
#    output:
#        "03.Alignment_hisat2/mapping.summary"
#    shell:
#        """
#	echo test
#        """
