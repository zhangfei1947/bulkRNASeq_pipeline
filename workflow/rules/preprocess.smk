localrules: make_links

rule make_links:
    input:
        r1 = lambda wildcards: config['samples'][wildcards.sample]['read1'],
        r2 = lambda wildcards: config['samples'][wildcards.sample]['read2']
    output:
        r1 = "01.Rawdata/{sample}/{sample}_R1.fq.gz",
        r2 = "01.Rawdata/{sample}/{sample}_R2.fq.gz"
    log:
        "logs/prep/{sample}.log"
    shell:
        """
        ln -sf {input.r1} {output.r1} 2> {log}
        ln -sf {input.r2} {output.r2} 2> {log}
        """