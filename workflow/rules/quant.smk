localrules: fc_summary

rule featurecounts:
    input:
        bams = expand("03.Alignment_hisat2/{sample}/{sample}.bam", sample=config['samples']),
        anno = config['genome']['annotation']
    output:
        "04.Quant_featureCounts/counts_raw.tsv"
    log:
        "logs/quant/featurecounts.log"
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 60,
        nodes = 1,
        ntasks_per_node = 4
    shell:
        """
        module load GCC/12.3.0 Subread/2.0.8

        featureCounts \
        -T {threads} \
        -a {input.anno} \
        -o {output} \
        -F GTF -t exon -g gene_id \
        -s 2 \
        -p \
        --countReadPairs \
        -O -M --fraction \
        -P -B \
        -d 40 \
        -D 800 \
        {input.bams} > {log} 2>&1
        """

rule fc_summary:
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