import pandas as pd

samples_path = '/juno/work/shah/users/chois7/tickets/crcorganoid98/136PC2/rnaseq/resources/short_samples.txt'
samples = [s.strip() for s in open(samples_path).readlines()]

for ix, sample in enumerate(samples):
    cnt_path = f'/juno/work/shah/users/chois7/tickets/crcorganoid98/136PC2/rnaseq/results/featurecounts/{sample}.counts.txt'
    cnt = pd.read_table(cnt_path, comment='#')
    df = cnt.iloc[:, [0, 6]]
    df.columns = ['gene', sample]
    df.set_index('gene', inplace=True)
    if ix == 0:
        merged = df.copy()
    if ix >= 1:
        merged = merged.join(df)

groups = ['ClbPpos', 'dClbP', 'W0']
for i in range(len(groups)):
    gi = groups[i]
    for j in range(i+1, len(groups)):
        gj = groups[j]
        out_table_path = f'/juno/work/shah/users/chois7/tickets/crcorganoid98/136PC2/rnaseq/results/deseq/{gi}_vs_{gj}.table.tsv'
        df = merged.loc[:, (merged.columns.str.count(gi) > 0) | (merged.columns.str.count(gj) > 0)]
        cols = df.columns
        cols = [c.replace('136PC2_', '') for c in cols]
        df.columns = cols
        df.to_csv(out_table_path, sep='\t', index=True)

        out_coldata_path = f'/juno/work/shah/users/chois7/tickets/crcorganoid98/136PC2/rnaseq/results/deseq/{gi}_vs_{gj}.coldata.tsv'
        keys = df.columns.tolist()
        values = ['Case' if s.count(gj) == 0 else 'Control' for s in keys]
        coldf = pd.DataFrame(columns=['gene', 'condition'])
        for k, v in zip(keys, values):
            coldf.loc[coldf.shape[0]] = [k.replace('136PC2_', ''), v]
        coldf.set_index('gene', inplace=True)
        coldf.to_csv(out_coldata_path, sep='\t', index=True)
    # break
