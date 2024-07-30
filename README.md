# 2024-nsurp

Automate NSURP commands to check results

## Setup:

First, `git clone` this repository and cd into the folder.

Download the files if necessary or softlink them into this directory

example soft link:
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



Install mamba environment (contains snakemake)
```
mamba env create -f environment.yml
mamba activate 2024-nsurp
```

## Run:

If on farm, open a tmux and then start an `srun` session with 10 CPU and 50G RAM
```
srun -p med2 --time=1-00:00:00 -c 10 --mem 50GB --pty bash
```

Do a snakemake dryrun:
```
snakemake -n
```

Run the workflow:
```
snakemake -c 4
```
