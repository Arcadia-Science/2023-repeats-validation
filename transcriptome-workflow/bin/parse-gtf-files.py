import argparse
import csv

def parse_tsv(input_file, output_file):
    seen_rows = set()  # To keep track of unique rows

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        for row in reader:
            # Check if the row contains 'CDS' in the third column
            if len(row) >= 3 and row[2] == 'CDS':
                # Parse the 'gene_id' and 'protein_id' values
                gene_id = None
                protein_id = None

                for field in row[8].split('; '):
                    if field.startswith('gene_id'):
                        gene_id = field.split('"')[1]
                    elif field.startswith('protein_id'):
                        protein_id = field.split('"')[1]

                # Check if this row is unique; if not, skip it
                if (gene_id, protein_id) not in seen_rows:
                    seen_rows.add((gene_id, protein_id))
                    writer.writerow([gene_id, protein_id])

def main():
    parser = argparse.ArgumentParser(description='Parse and filter a TSV file.')
    parser.add_argument('input_file', help='Input TSV file')
    parser.add_argument('output_file', help='Output TSV file')

    args = parser.parse_args()

    parse_tsv(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
