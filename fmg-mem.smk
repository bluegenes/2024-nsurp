import os
import csv
import math
import pandas as pd

# define some variables
# samples_file = "ralstonia-wort-sigs.txt"
# basename = "ralstonia-wort"

samples_file =  "ralstonia-branchwater-0.2.wort-sigs.txt"
basename = "branchwater-c0.2"

# samples_file = "ralstonia-wort-sigs.h1000.txt" # first 1000
# samples_file = "ralstonia-wort-sigs.h10.txt"
database_basename = "ralstonia32"
ralstonia_database = "databases/ralstonia32.zip"
database_tax = "databases/ralstonia32.lin-taxonomy.csv"
database_lingroups = "databases/ralstonia32.lingroups.csv"

out_dir = "output.ralstonia_sra"
logs_dir = f"{out_dir}/logs"

# now, get the sample names from the signature file paths
SAMPLES = [os.path.splitext(os.path.basename(path))[0] for path in open(samples_file)]
KSIZE = 21
SCALED = 1000

rule all:
    input:
        expand(f"{out_dir}/gather-k{{ksize}}-sc{{scaled}}/{{sample}}.gather.csv", sample=SAMPLES, ksize=KSIZE, scaled=SCALED),
        expand(f"{out_dir}/tax-k{{ksize}}-sc{{scaled}}/{{sample}}-x-{database_basename}.lingroup.tsv", sample=SAMPLES, ksize=KSIZE, scaled=SCALED),
        expand(f"{out_dir}/tax-k{{ksize}}-sc{{scaled}}/{{sample}}-x-{database_basename}.summarized.csv", sample=SAMPLES, ksize=KSIZE, scaled=SCALED),
        expand(f"{out_dir}/{basename}-x-{database_basename}.k{{ksize}}-sc{{scaled}}.lingroup.tsv", ksize=KSIZE, scaled=SCALED),
        expand(f"{out_dir}/{basename}-x-{database_basename}.k{{ksize}}-sc{{scaled}}.best-phylogroup.tsv", ksize=KSIZE, scaled=SCALED),
        expand(f"{out_dir}/{basename}-x-{database_basename}.k{{ksize}}-sc{{scaled}}.all-phylogroup.tsv", ksize=KSIZE, scaled=SCALED),
        expand(f"{out_dir}/{basename}-x-{database_basename}.k{{ksize}}-sc{{scaled}}.phylogroup-compare.best-only.tsv", ksize=KSIZE, scaled=SCALED),
        expand(f"{out_dir}/{basename}-x-{database_basename}.k{{ksize}}-sc{{scaled}}.phylogroup-compare.all.tsv", ksize=KSIZE, scaled=SCALED),
        # expand(f"{out_dir}/gather/{{sample}}.gather.csv", sample=SAMPLES),
        # expand(f"{out_dir}/tax/{{sample}}-x-ralstonia32.lingroup.tsv", sample=SAMPLES),
        # expand(f"{out_dir}/tax/{{sample}}-x-ralstonia32.summarized.csv", sample=SAMPLES),
        # f"{out_dir}/{basename}-x-{database_basename}.lingroup.tsv",
        # f"{out_dir}/{basename}-x-{database_basename}.best-phylogroup.tsv",


# for mgx, in-memory fastmultigather takes ~80G memory for ~ 1000 files on 50 threads.
# tested 50 threads for 2000 files --> not much more memory (~3+ hrs)
# to run all: 50 threads, 150G memory --> ~1day18hrs
rule branchwater_fastmultigather:
    input:
        sample_file = samples_file,
        database = ralstonia_database,
    output:
        gather=expand(f"{out_dir}/gather-k{{ksize}}-sc{{scaled}}/{{sample}}.gather.csv", sample=SAMPLES, ksize=KSIZE, scaled=SCALED),
        # missing=expand(f"{out_dir}/gather-k{{ksize}}-sc{{scaled}}/gather.missing", ksize=KSIZE, scaled=SCALED),
    conda: "branchwater.yml"
    params:
        outd=f"{out_dir}/gather-k{KSIZE}-sc{SCALED}",
        ksize=KSIZE,
        scaled=SCALED,
    threads: 100
    log: f"{logs_dir}/fastmultigather.k{KSIZE}-sc{SCALED}.log"
    benchmark: f"{logs_dir}/fastmultigather.k{KSIZE}-sc{SCALED}.benchmark"
    shell:
        """
        sourmash scripts fastmultigather -k {params.ksize} \
                         --scaled {params.scaled} {input.sample_file} \
                         {input.database} >> fmg.missing 2> {log}
        mv *gather.csv *prefetch.csv {params.outd}
        """
        # ah - issue - if no results, gather file is not created. Need to touch it.

rule sourmash_taxonomy:
    input:
        gather_csv =  ancient(f"{out_dir}/gather-k{{ksize}}-sc{{scaled}}/{{sample}}.gather.csv"),
        taxonomy_csv = database_tax,
        lingroups_csv = database_lingroups
    output:
        out_lingroup=f"{out_dir}/tax-k{{ksize}}-sc{{scaled}}/{{sample}}-x-{database_basename}.lingroup.tsv",
        out_summary=f"{out_dir}/tax-k{{ksize}}-sc{{scaled}}/{{sample}}-x-{database_basename}.summarized.csv",
    params:
        output_base = f"{out_dir}/tax-k{{ksize}}-sc{{scaled}}/{{sample}}-x-{database_basename}"
    log: f"{logs_dir}/tax-k{{ksize}}-sc{{scaled}}/{{sample}}-x-{database_basename}.log"
    benchmark: f"{logs_dir}/tax-k{{ksize}}-sc{{scaled}}/{{sample}}-x-{database_basename}.benchmark"
    conda: "branchwater.yml"
    threads: 1
    shell:
        """
        sourmash tax metagenome --gather-csv {input.gather_csv} --lins \
                                     --taxonomy-csv {input.taxonomy_csv} \
                                     --lingroup {input.lingroups_csv} \
                                     -o {params.output_base} -F csv_summary 2> {log}
        """

# write all lingroup files to a text file that can be used to aggregate
rule write_lingroup_files:
    input:
        lingroup_files = expand(f"{out_dir}/tax-k{{ksize}}-sc{{scaled}}/{{sample}}-x-{database_basename}.lingroup.tsv", sample=SAMPLES, ksize=KSIZE, scaled=SCALED)
    output: f"{out_dir}/{basename}-x-{database_basename}.k{{ksize}}-sc{{scaled}}.lingroup_files.txt"
    run:
        with open(output[0], 'w') as f:
            for inF in input.lingroup_files:
                f.write(str(inF) + '\n')


# aggregate the lingroup files into a single file
rule aggregate_lingroup_taxonomy:
    input:
        lingroup_files = f"{out_dir}/{basename}-x-{database_basename}.k{{ksize}}-sc{{scaled}}.lingroup_files.txt"
    output: f"{out_dir}/{basename}-x-{database_basename}.k{{ksize}}-sc{{scaled}}.lingroup.tsv"
    threads: 1
    shell:
        """
        python aggregate-lingroups.py --output {output} --lingroups {input.lingroup_files}
        """


# now, extract the best phylogroup for each sample and merge with metadata
rule extract_phylogroup_and_merge_branchwater_metadata:
    input:
        agg_lingroups =f"{out_dir}/{basename}-x-{database_basename}.k{{ksize}}-sc{{scaled}}.lingroup.tsv",
        metadata = "ralstonia-branchwater.csv",
    output: 
        best=f"{out_dir}/{basename}-x-{database_basename}.k{{ksize}}-sc{{scaled}}.best-phylogroup.tsv",
        good=f"{out_dir}/{basename}-x-{database_basename}.k{{ksize}}-sc{{scaled}}.good-phylogroup.tsv",
        all=f"{out_dir}/{basename}-x-{database_basename}.k{{ksize}}-sc{{scaled}}.all-phylogroup.tsv",
    threads: 1
    shell:
        """
        python extract-best-match.py --lingroup-csv {input.agg_lingroups} \
                                     --all {output.all} \
                                     --best {output.best} \
                                     --good {output.good} \
                                     --metadata {input.metadata}
        """

rule compare_best_phylogroup:
    input:
        best_phylogroup = f"{out_dir}/{basename}-x-{database_basename}.k{{ksize}}-sc{{scaled}}.best-phylogroup.tsv",
        expected_phylogroups = "expected-phylogroup.csv",
    output: f"{out_dir}/{basename}-x-{database_basename}.k{{ksize}}-sc{{scaled}}.phylogroup-compare.best-only.tsv",
    threads: 1
    shell:
        """
        python compare-phylogroup.py --input {input.best_phylogroup} \
                                     --expected {input.expected_phylogroups} \
                                     --output {output}
        """


rule compare_all:
    input:
        all_phylogroup = f"{out_dir}/{basename}-x-{database_basename}.k{{ksize}}-sc{{scaled}}.all-phylogroup.tsv",
        expected_phylogroups = "expected-phylogroup.csv",
    output: f"{out_dir}/{basename}-x-{database_basename}.k{{ksize}}-sc{{scaled}}.phylogroup-compare.all.tsv"
    threads: 1
    shell:
        """
        python compare-phylogroup.py --input {input.all_phylogroup} \
                                     --expected {input.expected_phylogroups} \
                                     --output {output}
        """