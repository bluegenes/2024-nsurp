import argparse
import csv
import os

def filter_sigs(args):
    # Read the acc column from the CSV file
    with open(args.branchwater_csv, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        acc_set = {row['acc'] for row in reader}

    # Filter the sig file based on the acc_set
    with open(args.wort_sigs, 'r') as sigfile, open(args.output, 'w') as outfile:
        for line in sigfile:
            if os.path.basename(line.strip()).replace('.sig', '') in acc_set:
                outfile.write(line)

def main():
    parser = argparse.ArgumentParser(description="Filter SIG file based on CSV 'acc' column")
    parser.add_argument("--wort-sigs", help="Path to the SIG file")
    parser.add_argument("--branchwater-csv", help="Path to the CSV file")
    parser.add_argument('-o', "--output", help="Path to the output file")
    
    args = parser.parse_args()
    
    filter_sigs(args)

if __name__ == "__main__":
    main()
