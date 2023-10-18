for sample in `cat resources/samples.txt`; do 
    cmd="bsub -J $sample -eo \
        logs/$sample.featureCounts.err \
        -oo logs/$sample.featureCounts.log \
        -M 16 \
        singularity run -B /juno -B /home ~/chois7/singularity/sif/featurecount.sif \
        featureCounts -p -O -T 4 -a reference/gencode.v31.annotation.gtf \
        -g gene_name -o results/featurecounts/${sample}.counts.txt \
        results/star/${sample}/${sample}.Aligned.sortedByCoord.out.bam"; 
    echo $cmd; 
done
