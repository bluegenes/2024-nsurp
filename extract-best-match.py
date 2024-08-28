import argparse
import pandas as pd

def extract_and_filter_phylogroup(lin):
    if lin == 'O_outgroup':
        return 'outgroup'
    parts = lin.split(';')
    # only keep phylogroup rows, ignoring more specific rows
    if parts[-1].startswith('B_'):
    # if len(parts) == 2:
        return parts[-1].replace('B_', '')
    else:
        return None


def main(args):
    df = pd.read_csv(args.lingroup_csv, sep='\t')
    # ensure 'percent_containment' is considered numeric
    df['percent_containment'] = pd.to_numeric(df['percent_containment'])

    # Extract best match per sample at the phylogroup level

    # create phylogroup column
    df['phylogroup'] = df['name'].apply(extract_and_filter_phylogroup)
    pg = df.dropna()

    # group rows by 'sample' and get the row with the maximum 'percent_containment' for each group
    best_matches = pg.loc[pg.groupby(['sample'])['percent_containment'].idxmax()]

    above_threshold = pg[pg['percent_containment'] >= args.containment_threshold]
    total_samples = df['sample'].nunique()
    samples_above_threshold = above_threshold['sample'].nunique()
    percentage_over_1 = (samples_above_threshold / total_samples) * 100
    print(f"Percentage of samples with at least one containment value over {args.containment_threshold}%: {percentage_over_1:.2f}%")

   # concatenate and drop duplicates
    # good = pd.concat([best_matches, above_threshold]).drop_duplicates().reset_index(drop=True) 
    if args.metadata is not None:
        branchwater_metadata = pd.read_csv(args.metadata)
        best_matches = best_matches.merge(branchwater_metadata, left_on='sample', right_on="acc", how='left')
        above_threshold = above_threshold.merge(branchwater_metadata, left_on='sample', right_on="acc", how='left')

    # rename + drop some columns
    best_matches.rename(columns={'percent_containment': 'phylogroup_containment',
                                 'containment': 'branchwater_containment',
                                 'organism': 'sra_organism'}, inplace=True)
                                #'num_bp_contained': 'num_bp_contained_phylogroup'
    best_matches.drop(columns=['name', 'lin', 'num_bp_contained', 'acc'], inplace=True)
    above_threshold.rename(columns={'percent_containment': 'phylogroup_containment',
                                    'containment': 'branchwater_containment',
                                    'organism': 'sra_organism'}, inplace=True)
                                    #'num_bp_contained': 'num_bp_contained_phylogroup'
    above_threshold.drop(columns=['name', 'lin', 'num_bp_contained', 'acc'], inplace=True)

    # save best to output file
    best_matches.to_csv(args.best, sep='\t', index=False)
    print(f"Best phylogroup match per sample saved to '{args.best.name}'")
    if args.good is not None:
        above_threshold.to_csv(args.good, sep='\t', index=False)
        print(f"Phylogroup matches above {args.containment_threshold}% containment saved to '{args.good.name}'")


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--lingroup-csv', type=argparse.FileType('r'), help='Aggregated lingroups input file', required=True)
    p.add_argument('-b', '--best', type=argparse.FileType('w'), help='Best phylogroup match', required=True)
    p.add_argument('--containment-threshold', type=float, default=0.3, help='Containment threshold for "good" matches')
    p.add_argument('-g', '--good', type=argparse.FileType('w'), help='Phylogroup matches above threshold')
    p.add_argument('--metadata', type=argparse.FileType('r'), help='branchwater metadata file')
    args = p.parse_args()
    main(args)