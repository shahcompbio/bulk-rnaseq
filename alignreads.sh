#! /bin/bash
#BSUB -J myjobMPI
#BSUB -n 80
#BSUB -R span[ptile=8]
#BSUB -R rusage[mem=6]
#BSUB -W 12:00
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -q cpuqueue
module load star/2.5.3a
for i in /data/ganesh/stoon/cat_fastqfiles/107C1_cc6_ClbPpos3_IGO*
do
    STAR --runMode alignReads \
        --outFilterMultimapNmax 3 \
        --runThreadN 80 \
        --genomeDir /data/ganesh/genome \
        --readFilesIn $i ${i%R1_001.fastq.gz}R2_001.fastq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix ${i%R1_001.fastq.gz}. \
        --outSAMtype BAM SortedByCoordinate \
    done
exit
