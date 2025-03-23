import os
import argparse
from Bio import SeqIO
from collections import defaultdict

def merge_faa_files(input_dirs, output_dir):
    """
    Merges sequences from matching .faa files in multiple directories and saves them to an output directory.

    Args:
        input_dirs (list of str): List of input directory paths.
        output_dir (str): Path to the output directory.
    """

    try:
        os.makedirs(output_dir, exist_ok=True)  # Create output directory if it doesn't exist

        file_dict = defaultdict(list)

        # Collect files from each input directory
        for input_dir in input_dirs:
            for filename in os.listdir(input_dir):
                if filename.endswith(".faa"):
                    file_path = os.path.join(input_dir, filename)
                    file_dict[filename].append(file_path)

        # Process matching files
        for filename, file_paths in file_dict.items():
            if len(file_paths) == len(input_dirs):  # Check if file exists in all input dirs
                output_file_path = os.path.join(output_dir, filename)
                merged_sequences = [] #Initialize merged_sequences here

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
                print(f"Warning: {filename} does not exist in all input directories.")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge sequences from matching .faa files in multiple directories.")
    parser.add_argument("input_dirs", nargs="+", help="List of input directory paths.")
    parser.add_argument("output_dir", help="Path to the output directory.")

    args = parser.parse_args()

    merge_faa_files(args.input_dirs, args.output_dir)