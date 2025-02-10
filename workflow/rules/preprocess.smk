localrules: make_links

rule make_links:
    input:
        r1 = lambda wildcards: os.path.join(
            config['samples'][wildcards.sample]['raw_dir'],
            f"{config['samples'][wildcards.sample]['raw_base']}_1.fq.gz"
        ),
        r2 = lambda wildcards: os.path.join(
            config['samples'][wildcards.sample]['raw_dir'],
            f"{config['samples'][wildcards.sample]['raw_base']}_2.fq.gz"
        )
    output:
        r1 = "01.Rawdata/{sample}_R1.fq.gz",
        r2 = "01.Rawdata/{sample}_R2.fq.gz"
    params:
        raw_base = lambda wildcards: config['samples'][wildcards.sample]['raw_base']
    shell:
        """
        ln -sf {input.r1} {output.r1}
        ln -sf {input.r2} {output.r2}
        """

