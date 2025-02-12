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
        nodes = 1,
        ntasks = 1,
        cpus_per_task= 6,
        mem_mb = 10000
    params:
        index = config['genome']['index'],
        splicesites = config['genome']['splicesites'],
        extra = config['align']
    shell:
        """
        module load GCC/13.2.0  OpenMPI/4.1.6 HISAT2/2.2.1 SAMtools/1.21

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

rule hisat2_summary:
    input:
        expand("03.Alignment_hisat2/{sample}/{sample}.summary", sample=config['samples'])
    output:
        "03.Alignment_hisat2/mapping.summary"
    shell:
        """
        # Define output file name
        output_file="{output[0]}"
        
        # Initialize output file and add header
        echo "Sample Name\tTotal Reads\tExact Match Rate\tMultiple Match Rate" > $output_file
        
        # Loop through all summary files
        for file in */*.summary; do
            # Extract sample name
            sample_name=$(basename "$file" | sed 's/\.summary//')
            
            # Extract total reads
            total_reads=$(grep -oP '^\d+ reads' "$file" | awk '{print $1}')
            
            # Extract exact match rate
            exact_match_rate=$(grep -oP '(?<=aligned concordantly exactly 1 time)\s+\d+(\.\d+)%' "$file")
            
            # Extract multiple match rate
            multiple_match_rate=$(grep -oP '(?<=aligned concordantly >1 times)\s+\d+(\.\d+)%' "$file")
            
            # Append extracted data to output file
            echo "$sample_name\t$total_reads\t$exact_match_rate\t$multiple_match_rate" >> $output_file
        done
        
        echo "Summarized files have been combined into $output_file"
        """

rule mapping_plot:
    input:
        "03.Alignment_hisat2/mapping.summary"
    output:
        "03.Alignment_hisat2/mappingrate.boxplot.png"
    shell:
        """
        module load GCC/12.2.0 OpenMPI/4.1.4 R/4.3.1
        export R_LIBS_USER="/scratch/group/lilab/software/R_library/4.3"

        Rscript -e '
        library(ggplot2)
        data <- read.table("{input}", header=TRUE, sep = "\t")
        # Create boxplot
        p <- ggplot(data, aes(x="", y=Exact.Match.Rate)) +
          geom_boxplot(fill = "skyblue") +
          theme_minimal() +
          labs(title = "Boxplot of Exact Match Rate",
               x = "",
               y = "Exact Match Rate (%)") +
          theme(axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                axis.text.x = element_text(angle = 45, hjust = 1))
        # Add text labels for the three lowest points
        outliers <- subset(data, Exact.Match.Rate %in% c(min(Exact.Match.Rate, na.rm = TRUE), 
                                                       min(Exact.Match.Rate, na.rm = TRUE),
                                                       min(Exact.Match.Rate, na.rm = TRUE)))
        p <- p + geom_text(aes(label=Sample.Name, y=Exact.Match.Rate), vjust=-0.5)  
        # Save the plot as PNG
        ggsave("{output}", plot=p, width=6, height=6, dpi=300)
        '
        """






