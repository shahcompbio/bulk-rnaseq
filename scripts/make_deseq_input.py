#!/usr/bin/env python3

import argparse
import os
import pandas as pd


def load_samples(samples_file):
    with open(samples_file) as f:
        return [line.strip() for line in f if line.strip()]


def merge_featurecounts(samples, counts_dir):
    merged = None
    for i, sample in enumerate(samples):
        cnt_path = os.path.join(counts_dir, f"{sample}.counts.txt")
        cnt = pd.read_table(cnt_path, comment='#')
        cnt.columns = ['gene', 'chrom', 'start', 'end', 'strand', 'length', sample]
        df = cnt[['gene', sample]]
        df.set_index('gene', inplace=True)
        merged = df if merged is None else merged.join(df)
    return merged


def write_deseq_inputs(merged, groups, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for i in range(len(groups)):
        gi = groups[i]
        for j in range(i + 1, len(groups)):
            gj = groups[j]
            df = merged.loc[:, (merged.columns.str.contains(gi)) | (merged.columns.str.contains(gj))]

            # Normalize column names
            cols = df.columns
            cols = [c.split('_')[1] for c in cols]
            df.columns = cols

            table_path = os.path.join(output_dir, f"{gi}_vs_{gj}.table.tsv")
            df.to_csv(table_path, sep='\t')

            keys = df.columns.tolist()
            values = ['Case' if gj not in s else 'Control' for s in keys]
            coldata = pd.DataFrame({'gene': [k for k in keys], 'condition': values})
            coldata.set_index('gene', inplace=True)

            coldata_path = os.path.join(output_dir, f"{gi}_vs_{gj}.coldata.tsv")
            coldata.to_csv(coldata_path, sep='\t')


def main():
    parser = argparse.ArgumentParser(description="Merge featureCounts results and prepare DESeq2 input files.")
    parser.add_argument('--samples_file', '-s', required=True, help="Path to the short_samples.txt file.")
    parser.add_argument('--counts_dir', '-i', required=True, help="Path to the featureCounts results directory.")
    parser.add_argument('--output_dir', '-o', required=True, help="Directory where DESeq2 input files will be written.")
    parser.add_argument('--groups', nargs='+', default=['ClbPpos', 'dClbP', 'W0'],
                        help="List of group prefixes to compare (default: ClbPpos dClbP W0).")

    args = parser.parse_args()

    samples = load_samples(args.samples_file)
    merged = merge_featurecounts(samples, args.counts_dir)
    write_deseq_inputs(merged, args.groups, args.output_dir)


if __name__ == '__main__':
    main()
