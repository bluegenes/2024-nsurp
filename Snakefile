####################################################
# sketch pre-downloaded metagenomes containing
# Ralstonia spp. + classify strain via LINs
####################################################
basename = "mgx6"
samples_csv=f"inputs/{basename}.samples.csv"
out_dir = "output.ralstonia"
logs_dir = f"{out_dir}/logs"

database_basename = "ralstonia32"
database_zip = "databases/ralstonia32.zip"
database_tax = "databases/ralstonia32.lin-taxonomy.csv"
database_lingroups = "databases/ralstonia32.lingroups.csv"

rule all:
    input:
        expand(f"{out_dir}/{basename}.zip", basename=basename),
        expand(f"{out_dir}/{basename}-x-{database_basename}.fmg.csv", basename=basename, database_basename=database_basename),
        expand(f"{out_dir}/{basename}-x-{database_basename}.lingroup.csv", basename=basename, database_basename=database_basename),


rule branchwater_manysketch:
    input: samples_csv
    output: "{basename}.zip"
    log: f"{logs_dir}/manysketch/{{basename}}.log"
    benchmark: f"{logs_dir}/manysketch/{{basename}}.benchmark"
    conda: "branchwater.yml"
    threads: 10
    shell:
        """
        sourmash scripts manysketch {input} -p k=21,k=31,k=51,scaled=1000,abund -o {output} 2> {log}
        """

rule branchwater_index_database:
    input: database_zip
    output: 
        db_current = f"{out_dir}/{{database_basename}}.rocksdb/CURRENT",
    log: f"{logs_dir}/index/{{database_basename}}.log"
    benchmark: f"{logs_dir}/index/{{database_basename}}.benchmark"
    params:
        db_dir = f"{out_dir}/{{database_basename}}.rocksdb",
    conda: "branchwater.yml"
    threads: 1
    shell:
        """
        sourmash scripts index {input} -m DNA -k 21 --scaled 1000 -o {params.db_dir} 2> {log}
        """

rule branchwater_fastmultigather:
    input: 
        samples = "{basename}.zip",
        database = f"{out_dir}/{{database_basename}}.rocksdb/CURRENT"
    output: f"{out_dir}/{{basename}}-x-{{database_basename}}.fmg.csv"
    params:
        db_dir = f"{out_dir}/{{database_basename}}.rocksdb"
    log: f"{logs_dir}/fastmultigather/{{basename}}-x-{{database_basename}}.fmg.log"
    benchmark: f"{logs_dir}/fastmultigather/{{basename}}-x-{{database_basename}}.fmg.benchmark"
    conda: "branchwater.yml"
    threads: 4
    shell:
        """
        sourmash scripts fastmultigather -k 21 --scaled 1000 {input.samples} {params.db_dir} -o {output} 2> {log}
        """
    
rule sourmash_taxonomy:
    input: 
        gather_csv = f"{out_dir}/{{basename}}-x-{{database_basename}}.fmg.csv",
        taxonomy_csv = database_tax,
        lingroups_csv = database_lingroups
    output: f"{out_dir}/{{basename}}-x-{{database_basename}}.lingroup.csv"
    params:
        output_base = f"{out_dir}/{{basename}}-x-{{database_basename}}"
    log: f"{logs_dir}/taxonomy/{{basename}}-x-{{database_basename}}.log"
    benchmark: f"{logs_dir}/taxonomy/{{basename}}-x-{{database_basename}}.benchmark"
    conda: "branchwater.yml"
    threads: 1
    shell:
        """
        sourmash taxonomy metagenome --gather-csv {input.gather_csv} \
                                     --taxonomy-csv {input.taxonomy_csv} \
                                     --lingroup {input.lingroups_csv} \
                                     -o {params.output_base} 2> {log}
        """
