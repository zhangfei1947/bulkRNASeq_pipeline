localrules: hisat2_summary, mapping_plot

rule hisat2_align:
    input:
        r1 = "02.QC_fastp/{sample}/{sample}_R1.clean.fastq.gz",
        r2 = "02.QC_fastp/{sample}/{sample}_R2.clean.fastq.gz"
    output:
        "03.Alignment_hisat2/{sample}/{sample}.bam",
        "03.Alignment_hisat2/{sample}/{sample}.summary",
        "03.Alignment_hisat2/{sample}/{sample}.bam.bai"
    log:
        "logs/align/{sample}.log"
    threads: 6
    resources:
        runtime = 120,
        cpus_per_task= 6,
        mem_mb = 10000
    params:
        index = config['genome']['index'],
        splicesites = config['genome']['splicesites']
    shell:
        """
hisat2 -x {params.index} --known-splicesite-infile {params.splicesites} \
    -p {threads} \
    --rna-strandness RF \
    --summary-file {output[1]} \
    --no-unal \
    -1 {input.r1} -2 {input.r2} 2> {log} |
samtools view -@ {threads} -Sb |
samtools sort -@ {threads} -o {output[0]}
samtools index {output[0]}
        """

rule hisat2_summary:
    input:
        expand("03.Alignment_hisat2/{sample}/{sample}.summary", sample=config["samples"].keys())
    output:
        "03.Alignment_hisat2/mapping.summary"
    container: None
    shell:
        """
echo "Sample Name\tTotal Reads\tExact Match Rate\tMultiple Match Rate" > {output}
for file in 03.Alignment_hisat2/*/*.summary; 
do
    sample_name=$(basename "$file" | sed 's/.summary//')   
    total_reads=$(grep -oP '^\d+ reads' "$file" | cut -d " " -f1)  
    exact_match_rate=$(grep "aligned concordantly exactly 1 time" "$file" | grep -oP "\d+\.\d+")
    multiple_match_rate=$(grep "aligned concordantly >1 times" "$file" | grep -oP "\d+\.\d+")
    echo "$sample_name\t$total_reads\t$exact_match_rate\t$multiple_match_rate" >> {output}
done
        """

rule mapping_plot:
    input:
        "03.Alignment_hisat2/mapping.summary"
    output:
        "03.Alignment_hisat2/mappingrate.boxplot.png"
    script:
        "../scripts/maprate.plot.R"






