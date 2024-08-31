import os
import argparse
import csv

def extract_sample_name(filename):
    fn = os.path.basename(filename)
    return fn.split('-x-')[0]

def main(args):

    # read in lingroup file list
    with args.lingroups as infile:
        # read in the list of lingroup files and strip
        lg_files = [line.strip() for line in infile]

    with args.output as outfile:
            writer = csv.writer(outfile, delimiter='\t')
            # make sure we write the header only once
            header_written = False
            
            # Iterate over files in the directory
            for inF in lg_files:
                sample_name = extract_sample_name(str(inF))
                # open the input file
                with open(inF, mode='r') as inF:
                    reader = csv.reader(inF, delimiter='\t')
                    # write the header only once
                    if not header_written:
                        header = next(reader)
                        header.insert(0, 'sample') # prepend sample_name
                        writer.writerow(header)
                        header_written = True
                    else:
                        next(reader)  # Skip the header if already written
                        
                    # Write the data rows
                    for row in reader:
                        row.insert(0, sample_name) # prepend sample
                        writer.writerow(row)

    print(f"Combined file saved to '{args.output.name}'")


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--lingroups', type=argparse.FileType('r'), help='file with list of lingroup files', required=True)
    p.add_argument('-o', '--output', type=argparse.FileType('w'), help='Aggregated lingroups output', required=True)
    args = p.parse_args()
    main(args)