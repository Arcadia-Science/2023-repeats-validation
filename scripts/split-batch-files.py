import os
import csv
import argparse

""""
This script takes in a directory of CSV files and assesses whether there are more than n rows (default 20), and if so creates batches of files with n number of rows
This was created because pysradb will error out after having to query SRA with many records in the same file
The batch CSVs maintain the same original number with a suffix for the number of file in that batch
It also skips empty files from the original record query since that happens occasionally 

"""

def split_csv(input_file, batch_size, output_directory):
    with open(input_file, 'r') as csv_file:
        reader = csv.reader(csv_file)
        try:
            header = next(reader)  # Attempt to read the header
        except StopIteration:
            print(f"Skipping empty file: {input_file}")
            return

        rows = list(reader)

    num_batches = (len(rows) + batch_size - 1) // batch_size

    for batch_num in range(num_batches):
        start_idx = batch_num * batch_size
        end_idx = (batch_num + 1) * batch_size
        batch_rows = rows[start_idx:end_idx]

        if not batch_rows:
            continue

        output_filename = f"{os.path.splitext(os.path.basename(input_file))[0]}_{batch_num + 1}.csv"
        output_file = os.path.join(output_directory, output_filename)

        with open(output_file, 'w', newline='') as batch_csv:
            writer = csv.writer(batch_csv)
            writer.writerow(header)
            writer.writerows(batch_rows)

        print(f"Created {output_file} with {len(batch_rows)} rows.")

def main():
    parser = argparse.ArgumentParser(description="Split CSV files into batches.")
    parser.add_argument("input_directory", help="Directory containing CSV files")
    parser.add_argument("--batch_size", type=int, default=20, help="Maximum number of rows per batch")
    parser.add_argument("--output_directory", default="output", help="Directory to store the output CSV files")

    args = parser.parse_args()

    input_directory = args.input_directory
    batch_size = args.batch_size
    output_directory = args.output_directory

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for filename in os.listdir(input_directory):
        if filename.endswith('.csv'):
            input_file = os.path.join(input_directory, filename)
            split_csv(input_file, batch_size, output_directory)

if __name__ == "__main__":
    main()
