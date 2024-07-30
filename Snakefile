####################################################
# sketch pre-downloaded metagenomes containing
# Ralstonia spp. + classify strain via LINs
####################################################
import os
import pandas as pd

basename = "mgx6"
samples_csv=f"inputs/{basename}.samples.csv"
out_dir = "output.ralstonia"
logs_dir = f"{out_dir}/logs"

# read metagenome info to get sample names
mgx_info = pd.read_csv(samples_csv)
SAMPLES = mgx_info["name"].tolist() 
print(SAMPLES)

database_basename = "ralstonia32"
database_zip = "databases/ralstonia32.zip"
database_tax = "databases/ralstonia32.lin-taxonomy.csv"
database_lingroups = "databases/ralstonia32.lingroups.csv"

rule all:
    input:
        expand(f"{out_dir}/{basename}.zip", basename=basename),
        expand(f"{out_dir}/{{sample}}.lingroup.tsv", sample=SAMPLES)


rule branchwater_manysketch:
    input: samples_csv
    output: f"{out_dir}/{{basename}}.zip"
    log: f"{logs_dir}/manysketch/{{basename}}.log"
    benchmark: f"{logs_dir}/manysketch/{{basename}}.benchmark"
    conda: "branchwater.yml"
    threads: 10
    shell:
        """
        sourmash scripts manysketch {input} -p k=21,k=31,k=51,scaled=1000,abund -o {output} 2> {log}
        """

rule branchwater_fastmultigather:
    input: 
        samples = f"{out_dir}/{basename}.zip",
        database = database_zip
    output: expand(f"{out_dir}/{{sample}}.gather.csv", sample=SAMPLES)
    log: f"{logs_dir}/fastmultigather/{basename}.fmg.log"
    benchmark: f"{logs_dir}/fastmultigather/{basename}.fmg.benchmark"
    conda: "branchwater.yml"
    threads: 10
    shell:
        """
        sourmash scripts fastmultigather -k 21 --scaled 1000 {input.samples} {input.database} 2> {log}
        """

rule sourmash_taxonomy:
    input: 
        gather_csv = f"{out_dir}/{{sample}}.gather.csv",
        taxonomy_csv = database_tax,
        lingroups_csv = database_lingroups
    output: f"{out_dir}/{{sample}}.lingroup.tsv"
    params:
        output_base = f"{out_dir}/{{sample}}-x-{database_basename}"
    log: f"{logs_dir}/taxonomy/{{sample}}-x-{database_basename}.log"
    benchmark: f"{logs_dir}/taxonomy/{{sample}}-x-{database_basename}.benchmark"
    conda: "branchwater.yml"
    threads: 1
    shell:
        """
        sourmash tax metagenome --gather-csv {input.gather_csv} --lins \
                                     --taxonomy-csv {input.taxonomy_csv} \
                                     --lingroup {input.lingroups_csv} \
                                     -o {params.output_base} 2> {log}
        """
