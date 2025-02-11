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
        mem_mb = 4000,
        runtime = 30,
        ntasks = 1，
        cpus_per_task = 2
    shell:
        """
        ml
        module load GCC/11.2.0 fastp/0.23.2
        ml
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
        alias jq="/scratch/group/lilab/software/jq-linux-i386"

        echo -e "Sample\tReads_brfore\tReads_after\tQ20_R1_before\tQ20_R2_before\tQ20_R1_after\tQ20_R2_after\tduplication_rate\tinsert_size" > {output[0]}
        for f in {input}; do
            sample=$(basename $f .json)
            total_reads_b=$(jq '.summary.before_filtering.total_reads' $f)
            total_reads_a=$(jq '.summary.after_filtering.total_reads' $f)
            q20b_r1=$(jq '.summary.before_filtering.q20_rate' $f)
            q20b_r2=$(jq '.summary.before_filtering.q20_rate' $f)
            q20a_r1=$(jq '.summary.after_filtering.q20_rate' $f)
            q20a_r2=$(jq '.summary.after_filtering.q20_rate' $f)
            dup_rate=$(jq '.summary.duplication.rate' $f)
            insert_size=$(jq '.summary.insert_size.peak' $f)
            echo -e "$sample\t$total_reads_b\t$total_reads_a\t$q20b_r1\t$q20b_r2\t$q20a_r1\t$q20a_r2\t$dup_rate\t$insert_size" >> {output[0]}
        done
        
        # Generate MultiQC report
        multiqc 02.QC_fastp/reports/ -o 02.QC_fastp/
        """

rule duprate_plot:
    input:
        "02.QC_fastp/QC_summary_table.tsv"
    output:
        "02.QC_fastp/duprate.boxplot.png.tsv"
    shell:
        """
        module load GCC/12.2.0 OpenMPI/4.1.4 R/4.3.1
        export R_LIBS_USER="/scratch/group/lilab/software/R_library/4.3"

        Rscript -e '
        # Load necessary libraries
        library(ggplot2)
        library(dplyr)
        
        # Read the data
        data <- read.delim("{input}", header = TRUE, sep = "\t")
        
        # Calculate the ranking of samples based on duplication rate and filter out the top 3
        data <- data %>%
          mutate(rank = rank(-duplication_rate)) %>%  # Rank samples by duplication rate in descending order
          filter(rank <= 3)  # Keep only the top 3 samples
        
        # Plot the boxplot
        ggplot(data, aes(x = "", y = duplication_rate, fill = Sample)) +  # Create a boxplot with empty x-axis labels
          geom_boxplot(width = 0.5) +  # Set the width of the boxes
          geom_text(aes(label = Sample), vjust = -1, position = position_jitter(width = 0.2, height = 0)) +  # Add sample names above the boxes
          theme_minimal() +  # Use a minimal theme
          labs(title = "Boxplot of Duplication Rate",  # Title of the plot
               x = "",  # Empty x-axis label
               y = "Duplication Rate") +  # Y-axis label
          theme(axis.text.x = element_blank(),  # Remove x-axis text
                axis.ticks.x = element_blank())  # Remove x-axis ticks
        
        # Save the plot as PNG file
        ggsave("{output}", width = 8, height = 6, dpi = 300)  # Save the plot with specified dimensions and resolution
        '

        """




