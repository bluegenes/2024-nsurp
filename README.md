# 2024-nsurp

Automate NSURP Ralstonia LIN-based taxonomic profiling

## Setup (on `farm`):

First, `git clone` this repository and cd into the folder:
```
git clone https://github.com/bluegenes/2024-nsurp
cd 2024-nsurp
```

Softlink the metagenome files into this directory
```
ln -s /group/ctbrowngrp4/mmerid/Ralstonia/rawdata/*fastq.gz ./inputs
```
> If you now do `ls inputs`, you should see the input fastq files.

Download the sourmash database. We have the associated taxonomy info in the `databases` folder already:
```
curl -JLO https://osf.io/download/wxtk3/
mv ralstonia.sc1000.zip databases/ralstonia32.zip

ls databases # look at the database files
```

Install mamba environment (contains snakemake, which we'll use to run)
```
mamba env create -f environment.yml
mamba activate 2024-nsurp
```

## Run:

Open a tmux and then start an `srun` session with 10 CPU and 50G RAM
```
srun -p med2 --time=1-00:00:00 -c 10 --mem 50GB --pty bash
```

Do a snakemake dryrun:
```
snakemake -n
```

Run the workflow:
```
snakemake -c 10
```


## To do plotting:

Create the plotting environment:
```
mamba env create -f plot-env.yml
```

```
mamba activate plot-env

# run the Rmarkdown file
Rscript -e "rmarkdown::render('plot-phylogroup-map.Rmd')"
```
