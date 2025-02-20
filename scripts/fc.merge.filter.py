import glob
import os
import pandas as pd
import sys

quant_dir = sys.argv[1]

def extract_sample_name(filepath):
    """Extracts the sample name from the file path."""
    parts = filepath.split("/")
    return parts[-2]

def merge_featurecounts(folder_path):
    """Merges featureCounts output files into a single table."""
    all_files = glob.glob(folder_path+ "/*/*.fC")
    if not all_files:
        print(f"No .fC files found in {folder_path}")
        return None

    all_data = {}

    for file in all_files:
        print(f"Processing file: {file}")  # Print current file being processed
        try:
            with open(file, 'r') as f:
                header_line = f.readline()
                second_line = f.readline()
                if not second_line.strip():
                    print(f"Warning: file {file} is empty or only has one line. Skipping.")
                    continue
                sample_name = extract_sample_name(file)
                if not sample_name:
                    print(f"Warning: Could not extract sample name from {file}. Skipping.")
                    continue

                df = pd.read_csv(file, sep='\t', skiprows=[0])
                all_data[sample_name] = df[['Geneid', df.columns[6]]]

        except (FileNotFoundError, pd.errors.EmptyDataError, IndexError) as e:
            print(f"Error processing file {file}: {e}")
            continue

    if not all_data:
        print("No valid data found to merge.")
        return None

    # Merge all dataframes
    merged_df = None
    for sample, df in all_data.items():
        df = df.rename(columns={df.columns[1]: sample})
        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on='Geneid', how='outer')

    return merged_df


merged_table = merge_featurecounts(quant_dir)

if merged_table is not None:
    #sort by sample name
    numeric_cols = merged_table.columns[1:]
    numeric_cols = sorted(numeric_cols)
    merged_table = merged_table[['Geneid']+numeric_cols]
    try:
        merged_table[numeric_cols] = merged_table[numeric_cols].astype(int)
    except ValueError as e:
        print(f"Warning: Non-numeric values found in count columns: {e}. Filtering may not work as expected.")
        #sys.exit(1)
    except TypeError as e:
        print(f"Warning: TypeError occurred. Check data types: {e}")
        #sys.exit(1)

    output_file = os.path.join(quant_dir, "merged_featurecounts.tsv")
    merged_table.to_csv(output_file, sep='\t', index=False)
    print(f"Merged table saved to {output_file}")

    # Filtering based on total read counts
    row_sums = merged_table[numeric_cols].sum(axis=1)
    num_samples = len(numeric_cols)
    ##filter average readcount > 1
    filtered_df = merged_table[row_sums >= num_samples]

    output_file_filtered = os.path.join(quant_dir, "merged_featurecounts_filtered.tsv") #more clear file name
    filtered_df.to_csv(output_file_filtered, sep='\t', index=False)
    print(f"Filtered table (total reads >= number of samples) saved to {output_file_filtered}")
