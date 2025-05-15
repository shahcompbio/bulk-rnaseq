import os
import subprocess
import pandas as pd

os.makedirs(config['log_dir'], exist_ok=True)
os.makedirs(config['tmp_dir'], exist_ok=True)

SAMPLES = [s.strip() for s in open(config['samplesfile']).readlines()]
GROUPS = ['ClbPpos', 'dClbP', 'W0']
PAIRS = ['ClbPpos_vs_dClbP', 'ClbPpos_vs_W0', 'dClbP_vs_W0']

rule all:
    input:
        expand('results/deseq/{pair}.DEG.csv', pair=PAIRS),

def _get_r1s(wildcards):
    df = pd.read_table(config['fastq_pathfile'])
    df = df[df['sample']==wildcards.sample]
    assert df.shape[0] == 1, df
    paths_str = df['R1_paths'].iloc[0]
    paths = paths_str.split(',')
    return paths

def _get_r2s(wildcards):
    df = pd.read_table(config['fastq_pathfile'])
    df = df[df['sample']==wildcards.sample]
    assert df.shape[0] == 1, df
    paths_str = df['R2_paths'].iloc[0]
    paths = paths_str.split(',')
    return paths

rule star_align:
    input:
        r1s = _get_r1s,
        r2s = _get_r2s,
    output:
        bam = 'results/star/{sample}.Aligned.sortedByCoord.out.bam'
        #bam = 'results/{sample}.bam'
        # bam = 'results/star/{sample}.bam',
    params:
        r1s = lambda w, input: ','.join(input.r1s),
        r2s = lambda w, input: ','.join(input.r2s),
        out_dir = 'results/star/{sample}',
        genome_dir = config['star_genome_dir'],
    threads: 16,
    resources:
        mem_mb = 61440,
    shell:
        'module load star/2.5.3a ; '
        'STAR --runMode alignReads '
        '--outFilterMultimapNmax 3 '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--readFilesIn {params.r1s} {params.r2s} '
        '--readFilesCommand zcat '
        '--outFileNamePrefix {params.out_dir}/{wildcards.sample}. '
        '--outSAMtype BAM SortedByCoordinate'

rule featurecounts:
    input:
        #bam = 'results/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam',
        bam = 'results/star/{sample}.bam',
    output:
        counts = 'results/featurecounts/{sample}.counts.txt',
    params:
        gtf = config['reference_gtf'],
    threads: 4,
    singularity: 'docker://alexgilgal/featurecount:latest',
    shell:
        'featureCounts -p -O -T {threads} '
        '-a {params.gtf} -g gene_name '
        '-o {output.counts} {input.bam}'

rule index_and_softlink:
    input:
        bam = 'results/star/{sample}.Aligned.sortedByCoord.out.bam',
        # bam = 'results/star/{sample}.bam',
    output:
        bam = 'results/star/{sample}.bam',
        bai = 'results/star/{sample}.bam.bai',
    params:
        bai = 'results/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai',
    shell:
        'samtools index {input.bam} && '
        'ln -s {input.bam} {output.bam} && '
        'ln -s {params.bai} {output.bai}'

rule make_deseq_input:
    input:
        expand('results/featurecounts/{sample}.counts.txt', sample=SAMPLES),
    output:
        table = expand('results/deseq/{pair}.table.tsv', pair=PAIRS),
        coldata = expand('results/deseq/{pair}.coldata.tsv', pair=PAIRS),
    params:
        samplesfile = config['samplesfile'],
        groups = ' '.join(GROUPS),
        counts_dir = 'results/featurecounts',
        output_dir = 'results/deseq',
    shell:
        'python scripts/make_deseq_input.py -s {params.samplesfile} '
        '-i {params.counts_dir} -o {params.output_dir} --groups {params.groups}'

rule run_deseq2:
    input:
        table = 'results/deseq_input/{pair}.table.tsv',
        coldata = 'results/deseq_input/{pair}.coldata.tsv',
    output:
        deg = 'results/deseq/{pair}.DEG.csv',
        norm = 'results/deseq/{pair}.norm.csv',
    params:
        output_dir = 'results/deseq'
    shell:
        'Rscript scripts/run_DESeq.R {input.table} {input.coldata} {params.output_dir}'
