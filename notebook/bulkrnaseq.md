---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.5
  kernelspec:
    display_name: Python 3.7.12 64-bit (conda)
    language: python
    name: python3712jvsc74a57bd0dc0fe05456373cce17991fe9e2e9264df4f4d1e972d77814a9700f21c9e7a8e2
---

# Settings

```python
import os
import re
import subprocess
import pandas as pd
import numpy as np
```

```python
import seaborn as sns
import matplotlib.pyplot as plt
```

# Data

```python
patient = '136PC2'
slug = 'RNA|TR|Bulk RNA'
meta_path = '/juno/work/shah/users/chois7/tickets/crcorganoid98/136PC2/rnaseq/resources/metadata.csv'
meta = pd.read_csv(meta_path)
```

```python
meta = meta[
    (meta['individual_identifier'] == patient) &
    (meta['technique_slug'] == slug) 
]
```

```python
meta.head()
```

```python
meta['sample_identifier']
```

```python
data_cols = ['sample', 'R1_paths', 'R2_paths']
df = pd.DataFrame(columns=data_cols)
dst_dir = '/juno/work/shah/users/chois7/tickets/crcorganoid98/136PC2/rnaseq/fastq'
for riw, row in meta.iterrows():
    sample = row['sample_identifier']
    dst_paths = {'R1':[], 'R2':[]}
    out_dir = f'{dst_dir}/{sample}'
    if not os.path.exists(f'{out_dir}'): 
        subprocess.run(f'mkdir -p {out_dir}', shell=True)
    fastq_paths = eval(row['raw_data'])
    for fastq_path in fastq_paths:
        src_dir, src_fname = os.path.split(fastq_path)
        lane, strand, _segment = re.search('(L\d\d\d)_(R\d)_(\d\d\d)', src_fname).groups()
        dst_fname = f'{sample}.fastq.gz'
        dst_fname = f'{sample}__{lane}_{strand}.fastq.gz'
        dst_path = f'{out_dir}/{dst_fname}'
        cmd = f'ln -s {fastq_path} {dst_path}'
        # subprocess.run(cmd, shell=True)
        dst_paths[strand].append(dst_path)
    field = [sample, ','.join(dst_paths['R1']), ','.join(dst_paths['R2'])]
    df.loc[df.shape[0]] = field
```

```python
df_path = '/juno/work/shah/users/chois7/tickets/crcorganoid98/136PC2/rnaseq/resources/fastq_paths.tsv'
df.to_csv(df_path, sep='\t', index=False)
```

```python
df
```

```python
eval(row['raw_data'])
```

# Preprocess


## prep data for DESeq2

```python
import pandas as pd
```

```python
root_dir = '/juno/work/shah/users/chois7/tickets/crcorganoid98/107C1/rnaseq'
```

```python
# samples_path = '/juno/work/shah/users/chois7/tickets/crcorganoid98/136PC2/rnaseq/resources/short_samples.txt'
samples_path = f'{root_dir}/resources/short_samples.txt'
samples = [s.strip() for s in open(samples_path).readlines()]
```

```python
for ix, sample in enumerate(samples):
    # cnt_path = f'/juno/work/shah/users/chois7/tickets/crcorganoid98/136PC2/rnaseq/results/featurecounts/{sample}.counts.txt'
    cnt_path = f'{root_dir}/results/featurecounts/{sample}.counts.txt'
    cnt = pd.read_table(cnt_path, comment='#')
    df = cnt.iloc[:, [0, 6]]
    df.columns = ['gene', sample]
    df.set_index('gene', inplace=True)
    if ix == 0:
        merged = df.copy()
    if ix >= 1:
        merged = merged.join(df)
```

```python
groups = ['ClbPpos', 'dClbP', 'W0']
for i in range(len(groups)):
    gi = groups[i]
    for j in range(i+1, len(groups)):
        gj = groups[j]
        out_table_path = f'{root_dir}/results/deseq/{gi}_vs_{gj}.table.tsv'
        df = merged.loc[:, (merged.columns.str.count(gi) > 0) | (merged.columns.str.count(gj) > 0)]
        cols = df.columns
        cols = [c.replace('136PC2_', '') for c in cols]
        cols = [c.replace('107C1_', '') for c in cols]
        df.columns = cols
        df.to_csv(out_table_path, sep='\t', index=True)

        out_coldata_path = f'{root_dir}/results/deseq/{gi}_vs_{gj}.coldata.tsv'
        keys = df.columns.tolist()
        values = ['Case' if s.count(gj) == 0 else 'Control' for s in keys]
        coldf = pd.DataFrame(columns=['gene', 'condition'])
        for k, v in zip(keys, values):
            coldf.loc[coldf.shape[0]] = [k.replace('136PC2_', ''), v]
        coldf.set_index('gene', inplace=True)
        coldf.to_csv(out_coldata_path, sep='\t', index=True)
    # break
```

# Expression comparison

```python
import gseapy as gp
```

```python
from collections import defaultdict
```

```python
raws = defaultdict(dict)
exprs = defaultdict(dict)
degs = defaultdict(dict)
gseas = defaultdict(dict)
```

## BULK

```python
modality = 'bulk'
```

### 107C1: ClbPpos vs dClbP

```python
experiment = '107C1'
wkdir = f'/juno/work/shah/users/chois7/tickets/crcorganoid98/{experiment}/rnaseq/results/deseq'

labels = ['ClbPpos', 'dClbP', 'W0']
case, ctrl = labels[0], labels[1]
deg_path = f'{wkdir}/{case}_vs_{ctrl}.DEG.csv'
df = pd.read_csv(deg_path)
columns = df.columns.tolist()
columns[0] = 'gene_name'
df.columns = columns

expr_path = f'{wkdir}/{case}_vs_{ctrl}.norm.csv'
expr = pd.read_csv(expr_path)
expr.columns = ['gene_name'] + expr.columns[1:].tolist()
expr.set_index('gene_name', inplace=True)
exprs[modality][experiment] = expr.copy()

raw_path = f'{wkdir}/{case}_vs_{ctrl}.table.tsv'
raw = pd.read_table(raw_path)
raw.columns = ['gene_name'] + raw.columns[1:].tolist()
raw.set_index('gene_name', inplace=True)
raws[modality][experiment] = raw.copy()
```

#### DEG

```python
padj_cutoff = 0.05
deg = df[df['padj'] < padj_cutoff]
deg.set_index('gene_name', inplace=True)
```

```python
# degs[f'{experiment} pks+/pks- ({modality})'] = deg.copy()
```

```python
degs[modality][experiment] = deg.copy()
```

```python
deg = degs['bulk']['107C1']
```

```python
deg
```

```python
degs['pseudobulk']['107C1'].head(20)
```

```python
sigdeg = deg[np.abs(deg['log2FoldChange']) >= np.log2(1.5)]
print(sigdeg.shape)
sigdeg
```

```python
for gene in sigdeg.index:
    print(gene)
```

#### GSEA

```python
classes = ['Case'] * 3 + ['Ctrl'] * 3
```

```python
%%time
modality, experiment = 'bulk', '107C1'
res = gp.gsea(
    data=exprs[modality][experiment],
    gene_sets="KEGG_2019_Human",
    # gene_sets = "MSigDB_Hallmark_2020",
    cls=classes,
) # row -> genes, column-> samples

gseas[modality][experiment] = res
```

```python
modality, experiment = 'bulk', '107C1'
res = gseas[modality][experiment]
sigres = res.res2d[res.res2d['FWER p-val'] < 0.1]
if sigres.shape[0] == 0:
    sigres = res.res2d.iloc[:10]
term = sigres.Term
# gp.gseaplot(res.ranking, term=term[i], **res.results[term[i]])
axs = res.plot(terms=term[:10])
```

### 136PC2: ClbPpos vs dClbP

```python
experiment = '136PC2'
wkdir = f'/juno/work/shah/users/chois7/tickets/crcorganoid98/{experiment}/rnaseq/results/deseq'

labels = ['ClbPpos', 'dClbP', 'W0']
case, ctrl = labels[0], labels[1]
deg_path = f'{wkdir}/{case}_vs_{ctrl}.DEG.csv'
df = pd.read_csv(deg_path)
columns = df.columns.tolist()
columns[0] = 'gene_name'
df.columns = columns

expr_path = f'{wkdir}/{case}_vs_{ctrl}.norm.csv'
expr = pd.read_csv(expr_path)
expr.columns = ['gene_name'] + expr.columns[1:].tolist()
expr.set_index('gene_name', inplace=True)
exprs[modality][experiment] = expr.copy()

raw_path = f'{wkdir}/{case}_vs_{ctrl}.table.tsv'
raw = pd.read_table(raw_path)
raw.columns = ['gene_name'] + raw.columns[1:].tolist()
raw.set_index('gene_name', inplace=True)
raws[modality][experiment] = raw.copy()
```

#### DEG

```python
padj_cutoff = 0.05
deg = df[df['padj'] < padj_cutoff]
deg.set_index('gene_name', inplace=True)
```

```python
# degs[f'{experiment} pks+/pks- ({modality})'] = deg.copy()
```

```python
degs[modality][experiment] = deg.copy()
```

```python
deg = degs['bulk']['136PC2']
```

```python
sigdeg = deg[np.abs(deg['log2FoldChange']) >= np.log2(1.5)]
print(sigdeg.shape)
sigdeg
```

```python
for gene in sigdeg.index:
    print(gene)
```

##### heatmap

```python
mtx_max = 1
mtx_min = -1
hashtag_values = {case: mtx_max, ctrl: mtx_min}
```

```python
i = 4
genes = ['hashtag'] + res.res2d.Lead_genes.iloc[i].split(";")
ax = gp.heatmap(
    df = df.loc[genes], #res.heatmat.loc[genes],
    z_score=None,
    title=res.res2d.Term.iloc[i],
    figsize=(6,0.5*df.loc[genes].index.shape[0]),
    cmap=plt.cm.viridis,
    xticklabels=False,
)
colors = {pair_classes[0]:'#1fa287', pair_classes[1]:'#471365'}
handle = [plt.plot([], [],
          color=colors[label], marker="s", ms=5, ls="")[0] for label in colors]
legend = ax.legend(handles=handle, labels=colors.keys(), title="hashtag", loc=(1.05, 0.75), frameon=False)
ax.add_artist(legend);
```

#### GSEA

```python
classes = ['Case'] * 3 + ['Ctrl'] * 3
```

```python
%%time
modality, experiment = 'bulk', '136PC2'
res = gp.gsea(
    data=exprs[modality][experiment],
    gene_sets="KEGG_2019_Human",
    # gene_sets = "MSigDB_Hallmark_2020",
    cls=classes,
) # row -> genes, column-> samples

gseas[modality][experiment] = res
```

```python
modality, experiment = 'bulk', '136PC2'
res = gseas[modality][experiment]
sigres = res.res2d[res.res2d['FWER p-val'] < 0.1]
print(sigres.shape)
if sigres.shape[0] < 10:
    sigres = res.res2d.iloc[:10]
term = sigres.Term
# gp.gseaplot(res.ranking, term=term[i], **res.results[term[i]])
axs = res.plot(terms=term[:10])
```

## Pseudo-BULK

```python
modality = 'pseudobulk'
```

### 107C1: ClbPpos vs dClbP

```python
experiment = '107C1'
wkdir = f'/juno/work/shah/users/chois7/tickets/crcorganoid98/{experiment}/pseudobulk/results/deseq'

labels = ['ClbPpos', 'dClbP', 'W0']
case, ctrl = labels[0], labels[1]
deg_path = f'{wkdir}/{case}_vs_{ctrl}.DEG.csv'
df = pd.read_csv(deg_path)
columns = df.columns.tolist()
columns[0] = 'gene_name'
df.columns = columns

expr_path = f'{wkdir}/{case}_vs_{ctrl}.norm.csv'
expr = pd.read_csv(expr_path)
expr.columns = ['gene_name'] + expr.columns[1:].tolist()
expr.set_index('gene_name', inplace=True)
exprs[modality][experiment] = expr.copy()

raw_path = f'{wkdir}/{case}_vs_{ctrl}.table.tsv'
raw = pd.read_table(raw_path)
raw.columns = ['gene_name'] + raw.columns[1:].tolist()
raw.set_index('gene_name', inplace=True)
raws[modality][experiment] = raw.copy()
```

#### DEG

```python
padj_cutoff = 0.05
deg = df[df['padj'] < padj_cutoff]
deg.set_index('gene_name', inplace=True)
```

```python
# degs[f'{experiment} pks+/pks- ({modality})'] = deg.copy()
```

```python
degs[modality][experiment] = deg.copy()
```

```python
deg = degs['pseudobulk']['107C1']
```

```python
sigdeg = deg[np.abs(deg['log2FoldChange']) >= np.log2(1.5)]
print(sigdeg.shape)
sigdeg
```

```python
for gene in sigdeg.index:
    print(gene)
```

#### GSEA

```python
classes = ['Case'] * 3 + ['Ctrl'] * 3
```

```python
%%time
modality, experiment = 'pseudobulk', '107C1'
res = gp.gsea(
    data=exprs[modality][experiment],
    gene_sets="KEGG_2019_Human",
    # gene_sets = "MSigDB_Hallmark_2020",
    cls=classes,
) # row -> genes, column-> samples

gseas[modality][experiment] = res
```

```python
sigres.shape
```

```python
modality, experiment = 'pseudobulk', '107C1'
res = gseas[modality][experiment]
sigres = res.res2d[res.res2d['FWER p-val'] < 0.1]
if sigres.shape[0] == 0:
    sigres = res.res2d.iloc[:10]
term = sigres.Term
# gp.gseaplot(res.ranking, term=term[i], **res.results[term[i]])
axs = res.plot(terms=term[:10])
```

### 136PC2: ClbPpos vs dClbP

```python
experiment = '136PC2'
wkdir = f'/juno/work/shah/users/chois7/tickets/crcorganoid98/{experiment}/pseudobulk/results/deseq'

labels = ['ClbPpos', 'dClbP', 'W0']
case, ctrl = labels[0], labels[1]
deg_path = f'{wkdir}/{case}_vs_{ctrl}.DEG.csv'
df = pd.read_csv(deg_path)
columns = df.columns.tolist()
columns[0] = 'gene_name'
df.columns = columns

expr_path = f'{wkdir}/{case}_vs_{ctrl}.norm.csv'
expr = pd.read_csv(expr_path)
expr.columns = ['gene_name'] + expr.columns[1:].tolist()
expr.set_index('gene_name', inplace=True)
exprs[modality][experiment] = expr.copy()

raw_path = f'{wkdir}/{case}_vs_{ctrl}.table.tsv'
raw = pd.read_table(raw_path)
raw.columns = ['gene_name'] + raw.columns[1:].tolist()
raw.set_index('gene_name', inplace=True)
raws[modality][experiment] = raw.copy()
```

#### DEG

```python
padj_cutoff = 0.05
deg = df[df['padj'] < padj_cutoff]
deg.set_index('gene_name', inplace=True)
```

```python
# degs[f'{experiment} pks+/pks- ({modality})'] = deg.copy()
```

```python
degs[modality][experiment] = deg.copy()
```

```python
deg
```

```python
deg = degs['pseudobulk']['136PC2']
```

```python
sigdeg = deg[np.abs(deg['log2FoldChange']) >= np.log2(1.5)]
print(sigdeg.shape)
sigdeg
```

```python
for gene in sigdeg.index:
    print(gene)
```

#### GSEA

```python
classes = ['Case'] * 3 + ['Ctrl'] * 3
```

```python
%%time
modality, experiment = 'pseudobulk', '136PC2'
res = gp.gsea(
    data=exprs[modality][experiment],
    gene_sets="KEGG_2019_Human",
    # gene_sets = "MSigDB_Hallmark_2020",
    cls=classes,
) # row -> genes, column-> samples

gseas[modality][experiment] = res
```

```python
modality, experiment = 'pseudobulk', '136PC2'
res = gseas[modality][experiment]
sigres = res.res2d[res.res2d['FWER p-val'] < 0.1]
print(sigres.shape)
if sigres.shape[0] <= 10:
    sigres = res.res2d.iloc[:10]
term = sigres.Term
# gp.gseaplot(res.ranking, term=term[i], **res.results[term[i]])
axs = res.plot(terms=term[:10])
```

```python
modality, experiment = 'pseudobulk', '136PC2'
# res = gseas[modality][experiment]
sigres = res.res2d[res.res2d['FWER p-val'] < 0.1]
term = sigres.Term
# gp.gseaplot(res.ranking, term=term[i], **res.results[term[i]])
axs = res.plot(terms=term[:10])
```

## expr volcanos

```python
exprs['bulk']['107C1'].head(2)
```

```python
exprs['bulk']['136PC2'].head(2)
```

```python
exprs['pseudobulk']['107C1'].head(2)
```

```python
exprs['pseudobulk']['136PC2'].head(2)
```

```python
degs['bulk']['107C1'].head(5)
```

```python
degs['bulk']['136PC2'].head(5)
```

```python
degs['pseudobulk']['107C1'].head(5)
```

```python
degs['pseudobulk']['136PC2'].head(5)
```

```python
deg = degs['pseudobulk']['136PC2']
```

```python
sigdeg = deg[np.abs(deg['log2FoldChange']) >= np.log2(1.5)]
print(sigdeg.shape)
sigdeg
```

```python
for gene in sigdeg.index:
    print(gene)
```

```python
np.log2(1.5)
```

## fc cmp betw bulk and pseudobulk


### 107C1

```python
cmp = degs['107C1 pks+/pks- (bulk)'].join(degs['107C1 pks+/pks- (pseudobulk)'], lsuffix='_bulk', rsuffix='_pseudobulk').dropna()
```

```python
cmp = degs['107C1 pks+/pks- (pseudobulk)'].join(degs['107C1 pks+/pks- (bulk)'], lsuffix='_pseudobulk', rsuffix='_bulk').dropna()
```

```python
cmp
```

```python
fig, ax = plt.subplots(figsize=(3,3))
plot_data = cmp.copy()
ax.scatter(plot_data['log2FoldChange_bulk'], plot_data['log2FoldChange_pseudobulk'])
for gene, row in cmp.iterrows():
    xcoord = row['log2FoldChange_bulk']
    ycoord = row['log2FoldChange_pseudobulk']
    text = gene
    ax.text(x=xcoord, y=ycoord, s=text)
```

```python
bulkdf
```

```python
bulkdf.loc['ZNF682']
```

```python
bulkdf.loc['TRPM6']
```

```python
cmpdf.loc['TRPM6']
```

```python
bpdf.loc['TRPM6']
```

```python
bpdf.loc['ZNF682']
```

```python
bpdf.loc['LCN2']
```

```python
cmpdf.loc['LINC01029']
```

### 136PC2

```python
cmp = degs['136PC2 pks+/pks- (bulk)'].join(degs['136PC2 pks+/pks- (pseudobulk)'], lsuffix='_bulk', rsuffix='_pseudobulk').dropna()
```

```python
cmp = degs['136PC2 pks+/pks- (pseudobulk)'].join(degs['136PC2 pks+/pks- (bulk)'], lsuffix='_pseudobulk', rsuffix='_bulk').dropna()
```

```python
cmp
```

```python
fig, ax = plt.subplots(figsize=(5,5))
plot_data = cmp.iloc[:5, :]
ax.scatter(plot_data['log2FoldChange_bulk'], plot_data['log2FoldChange_pseudobulk'])
```

```python
degs['107C1 pks+/pks- (pseudobulk)']
```

## expression cmp betw 107C1 and 136PC2


### bulk

```python
aliquot_ixs = [1, 2, 3]
pos_samples = [f'ClbPpos{i}' for i in aliquot_ixs]
neg_samples = [f'dClbP{i}' for i in aliquot_ixs]
```

```python
bulkdf = pd.DataFrame()
bulkdf['107C1 pks+'] = exprs['107C1'][pos_samples].mean(axis=1)
bulkdf['107C1 pks-'] = exprs['107C1'][neg_samples].mean(axis=1)
bulkdf['136PC2 pks+'] = exprs['136PC2'][pos_samples].mean(axis=1)
bulkdf['136PC2 pks-'] = exprs['136PC2'][neg_samples].mean(axis=1)
```

```python
color_map = {'107C1':'red', '136PC2':'orange'}
modality = 'bulk'
plot_data = bulkdf.copy()
n_row, n_col = 2, 3
fig, axes = plt.subplots(n_row, n_col, figsize=(9, 6), gridspec_kw={'hspace':0.5, 'wspace':0.5})
fig.suptitle('Bulk RNAseq normalized expression')
cnt = 0
for i in range(plot_data.shape[1]):
    for j in range(i+1, plot_data.shape[1]):
        x = cnt % n_row
        y = cnt // n_row
        cnt += 1
        ax = axes[x][y]
        columns = plot_data.columns
        col_i, col_j = columns[i], columns[j]
        ax.set(xscale="log", yscale="log")
        sns.scatterplot(data=plot_data, x=col_i, y=col_j, ax=ax,
                        alpha=1, s=1, edgecolor=None) 
        
        sample_i = col_i.split(' ')[0]
        sample_j = col_j.split(' ')[0]
        deg_i = degs[f'{sample_i} pks+/pks- ({modality})']
        if deg_i.shape[0] > 100:
            deg_i = deg_j.iloc[:100, :].copy()
        color = color_map[sample_i]
        genes_in_table = set(deg_i.index) & set(plot_data.dropna().index)
        expr_xs, expr_ys = [], []
        for gene in genes_in_table:
            expr_x = plot_data.loc[gene][col_i]
            expr_xs.append(expr_x)
            expr_y = plot_data.loc[gene][col_j]
            expr_ys.append(expr_y)
        ax.scatter(x=expr_xs, y=expr_ys, s=2, alpha=0.5, color=color)
        
        if sample_i != sample_j:
            deg_j = degs[f'{sample_j} pks+/pks- ({modality})']
            if deg_j.shape[0] > 100:
                deg_j = deg_j.iloc[:100, :].copy()
            color = color_map[sample_j]
            print(sample_i, sample_j, color)
            genes_in_table = set(deg_j.index) & set(plot_data.dropna().index)
            expr_xs, expr_ys = [], []
            for gene in genes_in_table:
                expr_x = plot_data.loc[gene][col_i]
                expr_xs.append(expr_x)
                expr_y = plot_data.loc[gene][col_j]
                expr_ys.append(expr_y)
            ax.scatter(x=expr_xs, y=expr_ys, s=2, alpha=0.5, color=color)

# plt.tight_layout()
```

### pseudobulk

```python
import seaborn as sns
import matplotlib.pyplot as plt
```

```python
aliquot_ixs = [1, 2, 3]
pos_samples = [f'ClbPpos{i}' for i in aliquot_ixs]
neg_samples = [f'dClbP{i}' for i in aliquot_ixs]
```

```python
cmpdf = pd.DataFrame()
cmpdf['107C1 pks+'] = exprs['107C1'][pos_samples].mean(axis=1)
cmpdf['107C1 pks-'] = exprs['107C1'][neg_samples].mean(axis=1)
cmpdf['136PC2 pks+'] = exprs['136PC2'][pos_samples].mean(axis=1)
cmpdf['136PC2 pks-'] = exprs['136PC2'][neg_samples].mean(axis=1)
```

```python
color_map = {'107C1':'pink', '136PC2':'purple'}
modality = 'pseudobulk'
plot_data = cmpdf.copy()
n_row, n_col = 2, 3
fig, axes = plt.subplots(n_row, n_col, figsize=(9, 6), gridspec_kw={'hspace':0.5, 'wspace':0.5})
fig.suptitle('scRNAseq normalized expression')
cnt = 0
for i in range(plot_data.shape[1]):
    for j in range(i+1, plot_data.shape[1]):
        x = cnt % n_row
        y = cnt // n_row
        cnt += 1
        ax = axes[x][y]
        columns = plot_data.columns
        col_i, col_j = columns[i], columns[j]
        ax.set(xscale="log", yscale="log")
        sns.scatterplot(data=plot_data, x=col_i, y=col_j, ax=ax,
                        alpha=1, s=1, edgecolor=None) 
        
        sample_i = col_i.split(' ')[0]
        sample_j = col_j.split(' ')[0]
        deg_i = degs[f'{sample_i} pks+/pks- ({modality})']
        if deg_i.shape[0] > 100:
            deg_i = deg_i.iloc[:100, :].copy()
        color = color_map[sample_i]
        genes_in_table = set(deg_i.index) & set(plot_data.dropna().index)
        expr_xs, expr_ys = [], []
        for gene in genes_in_table:
            expr_x = plot_data.loc[gene][col_i]
            expr_xs.append(expr_x)
            expr_y = plot_data.loc[gene][col_j]
            expr_ys.append(expr_y)
        ax.scatter(x=expr_xs, y=expr_ys, s=5, color=color, alpha=0.5)
        
        if sample_i != sample_j:
            deg_j = degs[f'{sample_j} pks+/pks- ({modality})']
            if deg_j.shape[0] > 100:
                deg_j = deg_j.iloc[:100, :].copy()
            color = color_map[sample_j]
            print(sample_i, sample_j, color)
            genes_in_table = set(deg_j.index) & set(plot_data.dropna().index)
            expr_xs, expr_ys = [], []
            for gene in genes_in_table:
                expr_x = plot_data.loc[gene][col_i]
                expr_xs.append(expr_x)
                expr_y = plot_data.loc[gene][col_j]
                expr_ys.append(expr_y)
            ax.scatter(x=expr_xs, y=expr_ys, s=5, color=color, alpha=0.5)
plt.tight_layout()
```

## joint: bulk vs pseudobulk

```python
bpdf = bulkdf.join(cmpdf, lsuffix=' (bulk)', rsuffix=' (pseudobulk)').dropna()
```

```python
conditions = ['107C1 pks+', '107C1 pks-', '136PC2 pks+', '136PC2 pks-']
```

```python
degs.keys()
```

```python
plot_data = bpdf.copy()
n_row, n_col = 2, 2
fig, axes = plt.subplots(n_row, n_col, figsize=(6, 6), gridspec_kw={'hspace':0.5, 'wspace':0.5})
fig.suptitle('sc-/bulk RNAseq normalized expression')
cnt = 0
for cnt, condition in enumerate(conditions):
    x = cnt % n_row
    y = cnt // n_row
    cnt += 1
    ax = axes[x][y]
    col_i, col_j = f'{condition} (bulk)', f'{condition} (pseudobulk)'
    ax.set(xscale="log", yscale="log")
    sns.scatterplot(data=plot_data, x=col_i, y=col_j, ax=ax,
                    alpha=1, s=1, edgecolor=None) 
    
    marker = 'o'
    modality = 'bulk'
    marker = marker_map[modality]
    color_map = {'107C1':'red', '136PC2':'orange'}
    sample_i = col_i.split(' ')[0]
    sample_j = col_j.split(' ')[0]
    deg_i = degs[f'{sample_i} pks+/pks- ({modality})']
    if deg_i.shape[0] > 100:
        deg_i = deg_i.iloc[:100, :].copy()
    color = color_map[sample_i]
    genes_in_table = set(deg_i.index) & set(plot_data.dropna().index)
    expr_xs, expr_ys = [], []
    for gene in genes_in_table:
        expr_x = plot_data.loc[gene][col_i]
        expr_xs.append(expr_x)
        expr_y = plot_data.loc[gene][col_j]
        expr_ys.append(expr_y)
    ax.scatter(x=expr_xs, y=expr_ys, s=5, color=color, alpha=0.5, marker=marker)

    if sample_i != sample_j:
        deg_j = degs[f'{sample_j} pks+/pks- ({modality})']
        if deg_j.shape[0] > 100:
            deg_j = deg_j.iloc[:100, :].copy()
        color = color_map[sample_j]
        print(sample_i, sample_j, color)
        genes_in_table = set(deg_j.index) & set(plot_data.dropna().index)
        expr_xs, expr_ys = [], []
        for gene in genes_in_table:
            expr_x = plot_data.loc[gene][col_i]
            expr_xs.append(expr_x)
            expr_y = plot_data.loc[gene][col_j]
            expr_ys.append(expr_y)
        ax.scatter(x=expr_xs, y=expr_ys, s=5, color=color, alpha=0.5, marker=marker)
        
    modality = 'pseudobulk'
    color_map = {'107C1':'pink', '136PC2':'purple'}
    sample_i = col_i.split(' ')[0]
    sample_j = col_j.split(' ')[0]
    deg_i = degs[f'{sample_i} pks+/pks- ({modality})']
    if deg_i.shape[0] > 100:
        deg_i = deg_i.iloc[:100, :].copy()
    color = color_map[sample_i]
    genes_in_table = set(deg_i.index) & set(plot_data.dropna().index)
    expr_xs, expr_ys = [], []
    for gene in genes_in_table:
        expr_x = plot_data.loc[gene][col_i]
        expr_xs.append(expr_x)
        expr_y = plot_data.loc[gene][col_j]
        expr_ys.append(expr_y)
    ax.scatter(x=expr_xs, y=expr_ys, s=5, color=color, alpha=0.5, marker=marker)

    if sample_i != sample_j:
        deg_j = degs[f'{sample_j} pks+/pks- ({modality})']
        if deg_j.shape[0] > 100:
            deg_j = deg_j.iloc[:100, :].copy()
        color = color_map[sample_j]
        print(sample_i, sample_j, color)
        genes_in_table = set(deg_j.index) & set(plot_data.dropna().index)
        expr_xs, expr_ys = [], []
        for gene in genes_in_table:
            expr_x = plot_data.loc[gene][col_i]
            expr_xs.append(expr_x)
            expr_y = plot_data.loc[gene][col_j]
            expr_ys.append(expr_y)
        ax.scatter(x=expr_xs, y=expr_ys, s=5, color=color, alpha=0.5, marker=marker)
plt.tight_layout()
```

```python
plot_data = bpdf.copy()
n_row, n_col = 4, 3
fig, axes = plt.subplots(n_row, n_col, figsize=(3*n_col, 3*n_row), gridspec_kw={'hspace':0.5, 'wspace':0.5})
fig.suptitle('sc-/bulk RNAseq normalized expression')
cnt = 0
for i in range(plot_data.shape[1]):
    for j in range(i+1, plot_data.shape[1]):
        col_i, col_j = columns[i], columns[j]
        if col_i.split(' (')[0] == col_j.split(' (')[0]: continue
        if col_i.split(' (')[1] == col_j.split(' (')[1]: continue
        x = cnt % n_row
        y = cnt // n_row
        cnt += 1
        ax = axes[x][y]
        columns = plot_data.columns
        ax.set(xscale="log", yscale="log")
        sns.scatterplot(data=plot_data, x=col_i, y=col_j, ax=ax,
                        alpha=1, s=1, edgecolor=None) 
plt.tight_layout()
```

# GSEA

```python
import blitzgsea as blitz
import pandas as pd

# list available gene set libraries in Enrichr
blitz.enrichr.print_libraries()

# use enrichr submodule to retrieve gene set library
# library = blitz.enrichr.get_library("GO_Biological_Process_2023")
library = blitz.enrichr.get_library("KEGG_2019_Human")
```

## example

```python
test['1'].hist(bins=30)
```

## BULK


### 107C1

```python
signature = df.iloc[:, [0, 2]].dropna()
signature.columns = ['0', '1']
```

```python
signature['1'].hist(bins=30)
```

```python
%%time
# run enrichment analysis
result = blitz.gsea(signature, library)
```

```python
result[result['fdr']< 0.05]
```

### 136PC2

```python
signature = df.iloc[:, [0, 2]].dropna()
signature.columns = ['0', '1']
```

```python
signature['1'].hist(bins=30)
```

```python
%%time
# run enrichment analysis
result = blitz.gsea(signature, library)
```

```python
df.head(2)
```

```python
result[result['fdr']< 0.05]
```

```python
result[result['fdr']< 0.05]
```

# OLD

```python
%%time
res = gp.gsea(
    data=expr,
    gene_sets="GO_Biological_Process_2021",
    # gene_sets = "MSigDB_Hallmark_2020",
    cls=classes,
) # row -> genes, column-> samples
```

```python
# df = pd.DataFrame(columns=['hashtag'])
# # df = df.join(res.heatmat)
# df['hashtag'] = list(pd.Series(res.heatmat.columns.map(pairdata.obs['Classification'])).replace(hashtag_values))
# df.index = res.heatmat.columns
# df = df.join(res.heatmat.T).T
```

```python
term = res.res2d.Term
# gp.gseaplot(res.ranking, term=term[i], **res.results[term[i]])
axs = res.plot(terms=term[:10])
```

```python
mtx_max = 1
mtx_min = -1
hashtag_values = {case: mtx_max, ctrl: mtx_min}
```

```python
i = 4
genes = ['hashtag'] + res.res2d.Lead_genes.iloc[i].split(";")
ax = gp.heatmap(
    df = df.loc[genes], #res.heatmat.loc[genes],
    z_score=None,
    title=res.res2d.Term.iloc[i],
    figsize=(6,0.5*df.loc[genes].index.shape[0]),
    cmap=plt.cm.viridis,
    xticklabels=False,
)
colors = {pair_classes[0]:'#1fa287', pair_classes[1]:'#471365'}
handle = [plt.plot([], [],
          color=colors[label], marker="s", ms=5, ls="")[0] for label in colors]
legend = ax.legend(handles=handle, labels=colors.keys(), title="hashtag", loc=(1.05, 0.75), frameon=False)
ax.add_artist(legend);
```
