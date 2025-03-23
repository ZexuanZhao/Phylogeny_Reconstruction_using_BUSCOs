import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse
import os

def extract_and_rename_sequences(full_table_path, translated_fasta_path, species_name, output_dir):
    """
    Extracts sequences based on 'Best gene' from 'full_table.tsv', renames headers, and saves each sequence to a separate FASTA file.

    Args:
        full_table_path (str): Path to the 'full_table.tsv' file.
        translated_fasta_path (str): Path to the 'translated_protein.fasta' file.
        species_name (str): The species name to use for renaming the headers.
        output_dir (str): Path to the output directory where FASTA files will be saved.
    """

    try:
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        df = pd.read_csv(full_table_path, sep='\t')
        single_genes = df[df['Status'] == 'Single']
        fasta_sequences = list(SeqIO.parse(translated_fasta_path, 'fasta'))

        for index, row in single_genes.iterrows(): #Iterate through rows to access both Gene and Best gene columns
            gene_name = row['Gene']
            best_gene = row['Best gene']

            found_sequence = None
            for record in fasta_sequences:
                if str(best_gene) in record.id:
                    found_sequence = record
                    break

            if found_sequence:
                new_record = SeqRecord(Seq(str(found_sequence.seq)), id=species_name, description="")
                output_fasta_path = os.path.join(output_dir, f"{gene_name}.faa") #Use Gene column for filename
                SeqIO.write([new_record], output_fasta_path, 'fasta')
                print(f"Sequence for '{gene_name}' written to {output_fasta_path}")
            else:
                print(f"Warning: Sequence for '{best_gene}' not found in the FASTA file.")

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract and rename sequences based on 'Best gene' values, saving each to a separate FASTA file.")
    parser.add_argument("full_table_path", help="Path to the full_table.tsv file.")
    parser.add_argument("translated_fasta_path", help="Path to the translated_protein.fasta file.")
    parser.add_argument("species_name", help="Species name to use for renaming headers.")
    parser.add_argument("output_dir", help="Path to the output directory for FASTA files.")

    args = parser.parse_args()

    extract_and_rename_sequences(args.full_table_path, args.translated_fasta_path, args.species_name, args.output_dir)