rule make_links:
    input:
        raw_dir = lambda wildcards: config['samples'][wildcards.sample]['raw_dir']
    output:
        "01.Rawdata/{sample}"
    shell:
        """
        ln -sf {input.raw_dir} {output}
        """

