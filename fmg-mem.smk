import os
import csv
import math

# define some variables
samples_file = "ralstonia-wort-sigs.txt"
basename = "ralstonia-wort"
# samples_file = "ralstonia-wort-sigs.h2000.txt"
database_basename = "ralstonia32"
ralstonia_database = "databases/ralstonia32.zip"
database_tax = "databases/ralstonia32.lin-taxonomy.csv"
database_lingroups = "databases/ralstonia32.lingroups.csv"

out_dir = "output.ralstonia_sra"
logs_dir = f"{out_dir}/logs"

# now, get the sample names from the signature file paths
SAMPLES = [os.path.splitext(os.path.basename(path))[0] for path in open(samples_file)]

rule all:
    input:
        # expand(f"{out_dir}/gather/{{sample}}.gather.csv", sample=SAMPLES),
        expand(f"{out_dir}/tax/{{sample}}-x-ralstonia32.lingroup.tsv", sample=SAMPLES),
        f"{out_dir}/{basename}-x-{database_basename}.lingroup.tsv"


# for mgx, in-memory fastmultigather takes ~80G memory for ~ 1000 files on 50 threads.
# tested 50 threads for 2000 files --> not much more memory (~3+ hrs)
# to run all: 50 threads, 150G memory --> ~1day18hrs
rule branchwater_fastmultigather:
    input:
        sample_file = samples_file,
        database = ralstonia_database,
    output:
        expand(f"{out_dir}/gather/{{sample}}.gather.csv", sample=SAMPLES),
    conda: "branchwater.yml"
    params:
        outd=f"{out_dir}/gather",
    threads: 100
    log: f"{logs_dir}/fastmultigather.log"
    benchmark: f"{logs_dir}/fastmultigather.benchmark"
    shell:
        """
        sourmash scripts fastmultigather -k 21 \
                         --scaled 1000 {input.sample_file} \
                         {input.database} 2> {log}
        mv *gather.csv *prefetch.csv {params.outd}
        """

rule sourmash_taxonomy:
    input: 
        gather_csv = f"{out_dir}/gather/{{sample}}.gather.csv",
        taxonomy_csv = database_tax,
        lingroups_csv = database_lingroups
    output: f"{out_dir}/tax/{{sample}}-x-{database_basename}.lingroup.tsv"
    params:
        output_base = f"{out_dir}/tax/{{sample}}-x-{database_basename}"
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

# aggregate the lingroup files into a single file
## for now, just aggregate all lines and pick later. Later, perhaps filter for best result??
rule aggregate_lingroup_taxonomy:
    input: 
        tax_files = expand(f"{out_dir}/tax/{{sample}}-x-{database_basename}.lingroup.tsv", sample=SAMPLES)
    output: f"{out_dir}/{basename}-x-{database_basename}.lingroup.tsv"
    log: f"{logs_dir}/taxonomy/{database_basename}.aggregate-lingroup.log"
    benchmark: f"{logs_dir}/taxonomy/{database_basename}.aggregate-lingroup.benchmark"
    threads: 1
    run: 
        def extract_sample_name(filename):
            fn = os.path.basename(filename)
            return fn.split('-x-')[0]

        with open(output_file, mode='w', newline='') as outfile:
            writer = csv.writer(outfile, delimiter='\t')
            # Write the header
            header_written = False
            
            # Iterate over files in the directory
            for filename in input.tax_files:
                sample_name = extract_sample_name(filename)
                # Open the input file
                with open(filename, mode='r') as inF:
                    reader = csv.reader(inF, delimiter='\t')
                    # Write the header only once
                    if not header_written:
                        header = next(reader)
                        # header.append('sample_name') 
                        header.insert(0, 'sample') # prepend sample_name
                        writer.writerow(header)
                        header_written = True
                    else:
                        next(reader)  # Skip the header
                        
                    # Write the data rows
                    for row in reader:
                        # row.append(sample_name)
                        row.insert(0, sample_name) # prepend sample
                        writer.writerow(row)

        print(f"Combined file saved to {output_file}")