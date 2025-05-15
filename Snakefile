import os
import subprocess
import pandas as pd

if not os.path.exists(config['log_dir']): subprocess.run(f'mkdir -p {config["log_dir"]}', shell=True)
if not os.path.exists(config['tmp_dir']): subprocess.run(f'mkdir -p {config["tmp_dir"]}', shell=True)


SAMPLES = [s.strip() for s in open(config['samplesfile_short']).readlines()]

rule all:
    input:
        #expand('results/star/{sample}.bam.bai', sample=SAMPLES),
        #expand('results/star/{sample}.bam.bai', sample=SAMPLES),
        expand('results/featurecounts/{sample}.counts.txt', sample=SAMPLES),
        #expand('results/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam', sample=SAMPLES),


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
        #bam = 'results/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam'
        #bam = 'results/{sample}.bam'
        bam = 'results/star/{sample}.bam',
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
    #singularity: 'docker://alexgilgal/featurecount:latest',
    singularity: '~/chois7/singularity/sif/featurecount.sif',
    shell:
        'featureCounts -p -O -T {threads} '
        '-a {params.gtf} -g gene_name '
        '-o {output.counts} {input.bam}'

rule index_and_softlink:
    input:
        #bam = 'results/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam',
        bam = 'results/star/{sample}.bam',
    output:
        bam = 'results/star/{sample}.bam',
        bai = 'results/star/{sample}.bam.bai',
    params:
        bai = 'results/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai',
    shell:
        'samtools index {input.bam}; '
        'cp {input.bam} {output.bam}; '
        'cp {params.bai} {output.bai}'

    
