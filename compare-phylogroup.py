import os
import pandas as pd

# function to check for phylogroup match
def check_phylogroup_match(row):
    if pd.notna(row['expected_phylogroup']):
        expected_groups = row['expected_phylogroup'].split(';')
        return row['phylogroup'] in expected_groups
    return False

def main(args):
    # Load the expected phylogroup file
    expected_df = pd.read_csv(args.expected_phylogroup)
    expected_df.rename(columns={'country': 'known_country'}, inplace=True)

    # Load the empirical results file
    detected_df = pd.read_csv(args.input, delimiter='\t')

    # Merge the two DataFrames on the 'geo_loc_name_country_calc' and 'country' columns
    merged_df = pd.merge(detected_df, expected_df, left_on='geo_loc_name_country_calc', right_on='known_country', how='left')

    # Create a new column 'expected_match' to check if the empirical phylogroup matches any of the expected phylogroups
    merged_df['expected_match'] = merged_df.apply(check_phylogroup_match, axis=1)

    # Save the result to a new CSV file
    merged_df.to_csv(args.output, index=False)

    print(f"Comparison result saved to '{args.output}'")
    # import pdb; pdb.set_trace()
    # to do - remove 'outgroup' from this count and also impose a threshold on the percent containment
    print(f"Total number of samples with detected phylogroup: {merged_df['phylogroup'].dropna().shape[0]}")

    # print the number of samples with a match compared to the number of samples with expected phylogroups
    num_matches = merged_df['expected_match'].sum()
    num_with_expected_phylogroups = merged_df['known_country'].dropna().shape[0]
    print (f"{num_with_expected_phylogroups} samples with expected phylogroups")
    print(f"phylogroups matching expectation: {num_matches}/{num_with_expected_phylogroups} ({num_matches/num_with_expected_phylogroups:.2%})")


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Compare detected phylogroups with expected phylogroups')
    parser.add_argument('--input', type=str, help='Path to the empirical phylogroup file')
    parser.add_argument('--expected-phylogroup', type=str, help='Path to the expected phylogroup file')
    parser.add_argument('--output', type=str, help='Path to save the comparison result')

    args = parser.parse_args()

    main(args)