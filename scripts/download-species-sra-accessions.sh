#!/bin/bash

# this documents the steps to download SRA record info for a given species or common species name
# the first part of this doesn't work well within a workflow because of how Entrez is installed and pulling down records
# is this pretty, not really. but does it work, yes.

# perform the esearch/efetch command with Entrez tools
# this will often fail for no good reason, so start back up in case the connection was interrupted etc.
for line in $(cat ../metadata/species_list.txt);
    do echo $line;
    esearch -db sra -query "$line" | efetch -format runinfo > $line-runinfo.csv;
done

# split records with more than 100 rows into batch files to pass to pysradb
# this has to be done because pysradb will error out after querying SRA for too long
python3 scripts/split-batch-files.py <INPUTDIR> --batch_size 100 --output_directory batchsrainfo

# get the detailed SRA info for each SRR accession in the runinfo with pysradb
# pass this to R filtering script to create one central metadata table
for file in *-runinfo.csv;
    do name=$(basename $file -runinfo.csv);
    pysradb metadata --saveto $name-srainfo.tsv --detailed $(cut -d "," -f1 $file | tail -n +2 | sort | uniq | xargs echo);
done
