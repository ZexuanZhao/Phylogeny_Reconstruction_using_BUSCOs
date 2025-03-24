import os
import argparse
from Bio import SeqIO
from collections import defaultdict
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

def merge_faa_files(input_dirs, output_dir, min_files, summary_table, num_threads):
    """
    Merges sequences from matching .faa files in multiple directories and saves them to an output directory.
    Creates a summary table of gene presence across directories.
    Runs the merging process in parallel using threads.

    Args:
        input_dirs (list of str): List of input directory paths.
        output_dir (str): Path to the output directory.
        min_files (int): Minimum number of files required for merging.
        summary_table (str): Path to the output summary table in TSV format.
        num_threads (int): Number of threads to use for parallel processing.
    """

    try:
        os.makedirs(output_dir, exist_ok=True)  # Create output directory if it doesn't exist

        file_dict = defaultdict(list)
        gene_presence = defaultdict(lambda: [0] * len(input_dirs))  # Initialize gene presence table

        # Collect files from each input directory and update gene presence
        for dir_index, input_dir in enumerate(input_dirs):
            for filename in os.listdir(input_dir):
                if filename.endswith(".faa"):
                    file_path = os.path.join(input_dir, filename)
                    file_dict[filename].append(file_path)
                    gene_presence[filename][dir_index] = 1  # Mark gene as present

        def process_file(filename, file_paths):
            if len(file_paths) >= min_files:  # Check if file exists in at least min_files input dirs
                output_file_path = os.path.join(output_dir, filename)
                merged_sequences = []  # Initialize merged_sequences here

                for file_path in file_paths:
                    try:
                        merged_sequences.extend(list(SeqIO.parse(file_path, "fasta")))
                    except FileNotFoundError:
                        print(f"Warning: File not found: {file_path}")
                        continue  # Skip to the next file path if one is not found.
                    except Exception as e:
                        print(f"An error occurred while reading {file_path}: {e}")
                        continue  # Skip to the next file path if an error occurs

                try:
                    SeqIO.write(merged_sequences, output_file_path, "fasta")
                    print(f"Merged sequences from {filename} to {output_file_path}")
                except Exception as e:
                    print(f"An error occurred while writing {output_file_path}: {e}")
            else:
                print(f"Warning: {filename} does not exist in at least {min_files} input directories.")

        # Process matching files in parallel
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            futures = [executor.submit(process_file, filename, file_paths) for filename, file_paths in file_dict.items()]
            for future in futures:
                future.result()  # Wait for all tasks to complete

        # Create and save summary table
        if summary_table:
            summary_data = {
                "gene": list(gene_presence.keys()),
                **{input_dirs[i]: [gene_presence[gene][i] for gene in gene_presence] for i in range(len(input_dirs))}
            }
            summary_df = pd.DataFrame(summary_data)
            try:
                summary_df.to_csv(summary_table, sep='\t', index=False)
                print(f"Summary table saved to {summary_table}")
            except Exception as e:
                print(f"An error occurred while saving summary table: {e}")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge sequences from matching .faa files in multiple directories.")
    parser.add_argument("input_dirs", nargs="+", help="List of input directory paths.")
    parser.add_argument("output_dir", help="Path to the output directory.")
    parser.add_argument("-n", "--min_files", type=int, default=1, help="Minimum number of files required for merging.")
    parser.add_argument("-s", "--summary_table", type=str, default=None, help="Path to the output summary table in TSV format.")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use for parallel processing.")

    args = parser.parse_args()

    merge_faa_files(args.input_dirs, args.output_dir, args.min_files, args.summary_table, args.threads)