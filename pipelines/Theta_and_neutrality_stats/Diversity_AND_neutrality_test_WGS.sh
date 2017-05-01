module load BCFtools/1.2-goolf-1.7.20
module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1

for i in *.bam;do samtools index $i;done
for i in *.fasta; do samtools faidx $i;done

ls Pj*.bam > Pj.bam.filelist
ls Pc*.bam > Pc.bam.filelist
ls Pm*.bam > Pm.bam.filelist

## WITHOUT ANCESTRAL STATES => OK
# Pj
~/utils/angsd/angsd -bam Pj.bam.filelist -doSaf 1 -anc GCA_001477535.1_Pneu_jiro_RU7_V2_genomic.fasta -GL 2 -P 5 -out Pj.outFold -fold 1
~/utils/angsd/misc/realSFS Pj.outFold.saf.idx -P 5  > Pj.outFold.saf.sfs
~/utils/angsd/angsd -bam Pj.bam.filelist -out Pj.outFold -doThetas 1 -doSaf 1 -pest Pj.outFold.saf.sfs -anc GCA_001477535.1_Pneu_jiro_RU7_V2_genomic.fasta -GL 2 -fold 1
~/utils/angsd/misc/thetaStat make_bed Pj.outFold.thetas.gz

#calculate Tajimas D
~/utils/angsd/misc/thetaStat do_stat Pj.outFold.thetas.gz -nChr 70 -win 50000 -step 10000  -outnames Pj.outFold.thetasWindow.gz

# Pc
~/utils/angsd/angsd -bam Pc.bam.filelist -doSaf 1 -anc GCA_001477545.1_Pneu_cari_B80_V3_genomic.fasta -GL 2 -P 5 -out Pc.outFold -fold 1
~/utils/angsd/misc/realSFS Pc.outFold.saf.idx -P 5 > Pc.outFold.saf.sfs
~/utils/angsd/angsd -bam Pc.bam.filelist -out Pc.outFold -doThetas 1 -doSaf 1 -pest Pc.outFold.saf.sfs -anc GCA_001477545.1_Pneu_cari_B80_V3_genomic.fasta -GL 2 -fold 1
~/utils/angsd/misc/thetaStat make_bed Pc.outFold.thetas.gz
~/utils/angsd/misc/thetaStat do_stat Pc.outFold.thetas.gz -nChr 62 -win 50000 -step 10000  -outnames Pc.outFold.thetasWindow.gz

# Pm
~/utils/angsd/angsd -bam Pm.bam.filelist -doSaf 1 -anc GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta -GL 2 -P 5 -out Pm.outFold -fold 1
~/utils/angsd/misc/realSFS Pm.outFold.saf.idx -P 5 > Pm.outFold.saf.sfs
~/utils/angsd/angsd -bam Pm.bam.filelist -out Pm.outFold -doThetas 1 -doSaf 1 -pest Pm.outFold.saf.sfs -anc GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta -GL 2 -fold 1
~/utils/angsd/misc/thetaStat make_bed Pm.outFold.thetas.gz
~/utils/angsd/misc/thetaStat do_stat Pm.outFold.thetas.gz -nChr 17 -win 50000 -step 10000  -outnames Pm.outFold.thetasWindow.gz