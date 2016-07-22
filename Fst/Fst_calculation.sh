#~/utils/last-744/src/lastdb -cR01 Pmdb GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta 

# SETTINGS
export LASTDIR=/nethome/cisseoh/utils/last-744
export READSDIR=/hpcdata/dcr_cc/Ousmane_Data/DATA/Pneumocystis_data/Reads
export REFGENOME=/hpcdata/dcr_cc/Ousmane_Data/DATA/Pneumocystis_data
export WRKDIR=/hpcdata/dcr_cc/Ousmane_Data/projects/Introgression/LAST_mapping/LAST_a15_b3_AND_LAST_Split

cd $WRKDIR

# LAST training 
$LASTDIR/src/lastdb -cR01 Pjtestdb GCA_001477535.1_Pneu_jiro_RU7_V2_genomic.fasta
$LASTDIR/scripts/last-train Pjdb GCA_001477545.1_Pneu_cari_B80_V3_genomic.fasta \
-P 4 -Q 0 > Pj_Pc_training.txt
# best scoring scheme:	lastal -j7 -S1 -Q0 -r5 -q5 -a15 -b3  # best model
$LASTDIR/scripts/last-train Pjdb GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta \
-P 8 -Q 0 > Pj_Pm_training.txt
# best scoring scheme:	-j7 -S1 -P8 -Q0 -r5 -q5 -a15 -b3 # best model

# DATA processing
# - P. jirovecii samples
# Pj RU817:
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 8 Pmdb \
$READSDIR/Pj_RU817_genome_reads/SRR1043750_1.fastq \
| $LASTDIR-744/src/last-split \
> Pj_RU817_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 8 Pmdb \
$READSDIR/Pj_RU817_genome_reads/SRR1043750_2.fastq \
| $LASTDIR-744/src/last-split \
> Pj_RU817_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pj_RU817_2_PmB123_1.maf > Pj_RU817_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pj_RU817_2_PmB123_2.maf > Pj_RU817_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools faidx GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pj_RU817_2_PmB123_1.sam > Pj_RU817_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pj_RU817_2_PmB123_2.sam > Pj_RU817_2_PmB123_2.Sam

samtools merge -uf Pj_RU817_2_PmB123.bam Pj_RU817_2_PmB123_1.Sam Pj_RU817_2_PmB123_2.Sam 
samtools sort Pj_RU817_2_PmB123.bam Pj_RU817_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
i=Pj_RU817_2_PmB123.sorted.bam \
O=Pj_RU817_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pj RU12
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 8 Pmdb \
$READSDIR/Pj_RU12_genome_reads/SRR1043747_1.fastq \
| $LASTDIR-744/src/last-split \
> Pj_RU12_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 8 Pmdb \
$READSDIR/Pj_RU12_genome_reads/SRR1043747_2.fastq \
> Pj_RU12_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pj_RU12_2_PmB123_1.maf > Pj_RU12_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pj_RU12_2_PmB123_2.maf > Pj_RU12_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pj_RU12_2_PmB123_1.sam > Pj_RU12_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pj_RU12_2_PmB123_2.sam > Pj_RU12_2_PmB123_2.Sam

samtools merge -uf Pj_RU12_2_PmB123.bam Pj_RU12_2_PmB123_1.Sam Pj_RU12_2_PmB123_2.Sam
samtools sort Pj_RU12_2_PmB123.bam Pj_RU12_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pj_RU12_2_PmB123.sorted.bam \
O=Pj_RU12_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pj Z
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 8 Pmdb \
$READSDIR/Pj_reads_Broad/PJ_Z_HighANDLow_merged_1.fq \
| $LASTDIR-744/src/last-split \
> Pj_Z_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 8 Pmdb \
$READSDIR/Pj_reads_Broad/PJ_Z_HighANDLow_merged_2.fq \
| $LASTDIR-744/src/last-split \
> Pj_Z_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pj_Z_2_PmB123_1.maf > Pj_Z_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pj_Z_2_PmB123_2.maf > Pj_Z_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pj_Z_2_PmB123_1.sam > Pj_Z_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pj_Z_2_PmB123_2.sam > Pj_Z_2_PmB123_2.Sam

samtools merge -uf Pj_Z_2_PmB123.bam Pj_Z_2_PmB123_1.Sam Pj_Z_2_PmB123_2.Sam
samtools sort Pj_Z_2_PmB123.bam Pj_Z_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pj_Z_2_PmB123.sorted.bam \
O=Pj_Z_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pj_W
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pj_reads_Broad/PJ_W_HighANDLow_merged_1.fq \
| $LASTDIR-744/src/last-split \
> Pj_W_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pj_reads_Broad/PJ_W_HighANDLow_merged_2.fq \
| $LASTDIR-744/src/last-split \
> Pj_W_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pj_W_2_PmB123_1.maf > Pj_W_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pj_W_2_PmB123_2.maf > Pj_W_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pj_W_2_PmB123_1.sam > Pj_W_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pj_W_2_PmB123_2.sam > Pj_W_2_PmB123_2.Sam

samtools merge -uf Pj_W_2_PmB123.bam Pj_W_2_PmB123_1.Sam Pj_W_2_PmB123_2.Sam
samtools sort Pj_W_2_PmB123.bam Pj_W_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pj_W_2_PmB123.sorted.bam \
O=Pj_W_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pj_SE8
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pj_SE8_genome_reads/Pj_SE8_Illumina_1.fastq \
| $LASTDIR-744/src/last-split \
> Pj_SE8_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pj_SE8_genome_reads/Pj_SE8_Illumina_2.fastq \
| $LASTDIR-744/src/last-split \
> Pj_SE8_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pj_SE8_2_PmB123_1.maf > Pj_SE8_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pj_SE8_2_PmB123_2.maf > Pj_SE8_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pj_SE8_2_PmB123_1.sam > Pj_SE8_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pj_SE8_2_PmB123_2.sam > Pj_SE8_2_PmB123_2.Sam

samtools merge -uf Pj_SE8_2_PmB123.bam Pj_SE8_2_PmB123_1.Sam Pj_SE8_2_PmB123_2.Sam
samtools sort Pj_SE8_2_PmB123.bam Pj_SE8_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pj_SE8_2_PmB123.sorted.bam \
O=Pj_SE8_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pj RN
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pj_RN110930_RNAseq/Pj_RN110930_RNAseq_merged_1.fastq \
| $LASTDIR-744/src/last-split \
> Pj_RN_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pj_RN110930_RNAseq/Pj_RN110930_RNAseq_merged_2.fastq \
| $LASTDIR-744/src/last-split \
> Pj_RN_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pj_RN_2_PmB123_1.maf > Pj_RN_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pj_RN_2_PmB123_2.maf > Pj_RN_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pj_RN_2_PmB123_1.sam > Pj_RN_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pj_RN_2_PmB123_2.sam > Pj_RN_2_PmB123_2.Sam
samtools merge -uf Pj_RN_2_PmB123.bam Pj_RN_2_PmB123_1.Sam Pj_RN_2_PmB123_2.Sam
samtools sort Pj_RN_2_PmB123.bam Pj_RN_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pj_RN_2_PmB123.sorted.bam \
O=Pj_RN_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pj RU7
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pj_RU7_genome_reads/SRR1043749_1.fastq \
| $LASTDIR-744/src/last-split \
> Pj_RU7_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pj_RU7_genome_reads/SRR1043749_2.fastq \
| $LASTDIR-744/src/last-split \
> Pj_RU7_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pj_RU7_2_PmB123_1.maf > Pj_RU7_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pj_RU7_2_PmB123_2.maf > Pj_RU7_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pj_RU7_2_PmB123_1.sam > Pj_RU7_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pj_RU7_2_PmB123_2.sam > Pj_RU7_2_PmB123_2.Sam
samtools merge -uf Pj_RU7_2_PmB123.bam Pj_RU7_2_PmB123_1.Sam Pj_RU7_2_PmB123_2.Sam
samtools sort Pj_RU7_2_PmB123.bam Pj_RU7_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pj_RU7_2_PmB123.sorted.bam \
O=Pj_RU7_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# ----- P. carinii
# Pc B80
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pc_B80_genome_reads/Pc_B80_genome_reads_merged_1.fastq \
| $LASTDIR-744/src/last-split \
> Pc_B80_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pc_B80_genome_reads/Pc_B80_genome_reads_merged_2.fastq \
| $LASTDIR-744/src/last-split \
> Pc_B80_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pc_B80_2_PmB123_1.maf > Pc_B80_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pc_B80_2_PmB123_2.maf > Pc_B80_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pc_B80_2_PmB123_1.sam > Pc_B80_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pc_B80_2_PmB123_2.sam > Pc_B80_2_PmB123_2.Sam
samtools merge -uf Pc_B80_2_PmB123.bam Pc_B80_2_PmB123_1.Sam Pc_B80_2_PmB123_2.Sam
samtools sort Pc_B80_2_PmB123.bam Pc_B80_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pc_B80_2_PmB123.sorted.bam \
O=Pc_B80_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pc B50
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pc_B50_genome_reads/SRR1043724_1.fastq \
| $LASTDIR-744/src/last-split \
> Pc_B50_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pc_B50_genome_reads/SRR1043724_2.fastq \
| $LASTDIR-744/src/last-split \
> Pc_B50_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pc_B50_2_PmB123_1.maf > Pc_B50_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pc_B50_2_PmB123_2.maf > Pc_B50_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pc_B50_2_PmB123_1.sam > Pc_B50_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pc_B50_2_PmB123_2.sam > Pc_B50_2_PmB123_2.Sam
samtools merge -uf Pc_B50_2_PmB123.bam Pc_B50_2_PmB123_1.Sam Pc_B50_2_PmB123_2.Sam
samtools sort Pc_B50_2_PmB123.bam Pc_B50_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pc_B50_2_PmB123.sorted.bam \
O=Pc_B50_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pc_B70
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pc_B70_genome_reads/SRR1043725_1.fastq \
| $LASTDIR-744/src/last-split \
> Pc_B70_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pc_B70_genome_reads/SRR1043725_2.fastq \
| $LASTDIR-744/src/last-split \
> Pc_B70_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pc_B70_2_PmB123_1.maf > Pc_B70_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pc_B70_2_PmB123_2.maf > Pc_B70_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pc_B70_2_PmB123_1.sam > Pc_B70_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pc_B70_2_PmB123_2.sam > Pc_B70_2_PmB123_2.Sam
samtools merge -uf Pc_B70_2_PmB123.bam Pc_B70_2_PmB123_1.Sam Pc_B70_2_PmB123_2.Sam
samtools sort Pc_B70_2_PmB123.bam Pc_B70_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pc_B70_2_PmB123.sorted.bam \
O=Pc_B70_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pc 70954
$LASTDIR-744/src/lastal -j7 -S1 -Q0 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pc_Broad_RNAseq/Pc_Solexa-70954_Left.fasta \
| $LASTDIR-744/src/last-split \
> Pc_Solexa-70954_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q0 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pc_Broad_RNAseq/Pc_Solexa-70954_Right.fasta \
| $LASTDIR-744/src/last-split \
> Pc_Solexa-70954_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pc_Solexa-70954_2_PmB123_1.maf > Pc_Solexa-70954_2_PmB123_1.sam 
$LASTDIR-744/scripts/maf-convert sam Pc_Solexa-70954_2_PmB123_2.maf > Pc_Solexa-70954_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pc_Solexa-70954_2_PmB123.sam > Pc_Solexa-70954_2_PmB123.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pc_Solexa-70954_2_PmB123_1.sam > Pc_Solexa-70954_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pc_Solexa-70954_2_PmB123_2.sam > Pc_Solexa-70954_2_PmB123_2.Sam

samtools merge -uf Pc_Solexa-70954_2_PmB123.bam Pc_Solexa-70954_2_PmB123_1.Sam Pc_Solexa-70954_2_PmB123_2.Sam
samtools sort Pc_Solexa-70954_2_PmB123.bam Pc_Solexa-70954_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pc_Solexa-70954_2_PmB123.sorted.bam \
O=Pc_Solexa-70954_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pc 70955
$LASTDIR-744/src/lastal -j7 -S1 -Q0 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pc_Broad_RNAseq/Pc_Solexa-70955_Left.fasta \
| $LASTDIR-744/src/last-split \
> Pc_Solexa-70955_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q0 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pc_Broad_RNAseq/Pc_Solexa-70955_Right.fasta \
| $LASTDIR-744/src/last-split \
> Pc_Solexa-70955_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pc_Solexa-70955_2_PmB123_1.maf > Pc_Solexa-70955_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pc_Solexa-70955_2_PmB123_2.maf > Pc_Solexa-70955_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pc_Solexa-70955_2_PmB123.sam > Pc_Solexa-70955_2_PmB123.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pc_Solexa-70955_2_PmB123_1.sam > Pc_Solexa-70955_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pc_Solexa-70955_2_PmB123_2.sam > Pc_Solexa-70955_2_PmB123_2.Sam

samtools merge -uf Pc_Solexa-70955_2_PmB123.bam Pc_Solexa-70955_2_PmB123_1.Sam Pc_Solexa-70955_2_PmB123_2.Sam
samtools sort Pc_Solexa-70955_2_PmB123.bam Pc_Solexa-70955_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pc_Solexa-70955_2_PmB123.sorted.bam \
O=Pc_Solexa-70955_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pc 70956
$LASTDIR-744/src/lastal -j7 -S1 -Q0 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pc_Broad_RNAseq/Pc_Solexa-70956_Left.fasta \
| $LASTDIR-744/src/last-split \
> Pc_Solexa-70956_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q0 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pc_Broad_RNAseq/Pc_Solexa-70956_Right.fasta \
| $LASTDIR-744/src/last-split \
> Pc_Solexa-70956_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pc_Solexa-70956_2_PmB123_1.maf > Pc_Solexa-70956_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pc_Solexa-70956_2_PmB123_2.maf > Pc_Solexa-70956_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pc_Solexa-70956_2_PmB123.sam > Pc_Solexa-70956_2_PmB123.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pc_Solexa-70956_2_PmB123_1.sam > Pc_Solexa-70956_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pc_Solexa-70956_2_PmB123_2.sam > Pc_Solexa-70956_2_PmB123_2.Sam

samtools merge -uf Pc_Solexa-70956_2_PmB123.bam Pc_Solexa-70956_2_PmB123_1.Sam Pc_Solexa-70956_2_PmB123_2.Sam
samtools sort Pc_Solexa-70956_2_PmB123.bam Pc_Solexa-70956_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pc_Solexa-70956_2_PmB123.sorted.bam \
O=Pc_Solexa-70956_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# ---- P. murina
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_B123_genome_reads/SRR770457_1.fastq \
| $LASTDIR-744/src/last-split \
> Pm_B123_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_B123_genome_reads/SRR770457_1.fastq \
| $LASTDIR-744/src/last-split \
> Pm_B123_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pm_B123_2_PmB123_1.maf > Pm_B123_2_PmB123_1.sam 
$LASTDIR-744/scripts/maf-convert sam Pm_B123_2_PmB123_2.maf > Pm_B123_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_B123_2_PmB123_1.sam > Pm_B123_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_B123_2_PmB123_2.sam > Pm_B123_2_PmB123_2.Sam

samtools merge -uf  Pm_B123_2_PmB123.bam Pm_B123_2_PmB123_1.Sam Pm_B123_2_PmB123_2.Sam
samtools sort Pm_B123_2_PmB123.bam Pm_B123_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pm_B123_2_PmB123.sorted.bam \
O=Pm_B123_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pm SAIC
$LASTDIR-744/src/lastal -Q0 -r5 -q5 -a35 -b5 -i8 -P 8 Pmdb \
$READSDIR/Pm_454_reads/P_murina_SAIC_all_runs_454Reads.fasta \
> Pm_saic_2_PmB123.maf
$LASTDIR-744/scripts/maf-convert sam Pm_saic_2_PmB123.maf > Pm_saic_2_PmB123.sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_saic_2_PmB123.sam > Pm_saic_2_PmB123.Sam
samtools view -Sb Pm_saic_2_PmB123.Sam > Pm_saic_2_PmB123.bam
samtools sort Pm_saic_2_PmB123.bam Pm_saic_2_PmB123.sorted
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pm_saic_2_PmB123.sorted.bam O=Pm_saic_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pm_Solexa-70950
$LASTDIR-744/src/lastal -j7 -S1 -Q0 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Broad_RNAseq/Pm_Solexa-70950_Left.fasta \
| $LASTDIR-744/src/last-split \
> Pm_Solexa-70950_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q0 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Broad_RNAseq/Pm_Solexa-70950_Right.fasta \
| $LASTDIR-744/src/last-split \
> Pm_Solexa-70950_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pm_Solexa-70950_2_PmB123_1.maf > Pm_Solexa-70950_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pm_Solexa-70950_2_PmB123_2.maf > Pm_Solexa-70950_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Solexa-70950_2_PmB123_1.sam > Pm_Solexa-70950_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Solexa-70950_2_PmB123_2.sam > Pm_Solexa-70950_2_PmB123_2.Sam

samtools merge -uf  Pm_Solexa-70950_2_PmB123.bam Pm_Solexa-70950_2_PmB123_1.Sam Pm_Solexa-70950_2_PmB123_2.Sam
samtools sort Pm_Solexa-70950_2_PmB123.bam Pm_Solexa-70950_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pm_Solexa-70950_2_PmB123.sorted.bam \
O=Pm_Solexa-70950_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pm_Solexa-70951
$LASTDIR-744/src/lastal -j7 -S1 -Q0 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Broad_RNAseq/Pm_Solexa-70951_Left.fasta \
| $LASTDIR-744/src/last-split \
> Pm_Solexa-70951_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q0 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Broad_RNAseq/Pm_Solexa-70951_Right.fasta \
| $LASTDIR-744/src/last-split \
> Pm_Solexa-70951_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pm_Solexa-70951_2_PmB123_1.maf > Pm_Solexa-70951_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pm_Solexa-70951_2_PmB123_2.maf > Pm_Solexa-70951_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Solexa-70951_2_PmB123_1.sam > Pm_Solexa-70951_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Solexa-70951_2_PmB123_2.sam > Pm_Solexa-70951_2_PmB123_2.Sam

samtools merge -uf  Pm_Solexa-70951_2_PmB123.bam Pm_Solexa-70951_2_PmB123_1.Sam Pm_Solexa-70951_2_PmB123_2.Sam
samtools sort Pm_Solexa-70951_2_PmB123.bam Pm_Solexa-70951_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pm_Solexa-70951_2_PmB123.sorted.bam \
O=Pm_Solexa-70951_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pm_Solexa-70952
$LASTDIR-744/src/lastal -j7 -S1 -Q0 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Broad_RNAseq/Pm_Solexa-70952_Left.fasta \
| $LASTDIR-744/src/last-split \
> Pm_Solexa-70952_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q0 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Broad_RNAseq/Pm_Solexa-70952_Right.fasta \
| $LASTDIR-744/src/last-split \
> Pm_Solexa-70952_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pm_Solexa-70952_2_PmB123_1.maf > Pm_Solexa-70952_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pm_Solexa-70952_2_PmB123_2.maf > Pm_Solexa-70952_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Solexa-70952_2_PmB123_1.sam > Pm_Solexa-70952_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Solexa-70952_2_PmB123_2.sam > Pm_Solexa-70952_2_PmB123_2.Sam

samtools merge -uf  Pm_Solexa-70952_2_PmB123.bam Pm_Solexa-70952_2_PmB123_1.Sam Pm_Solexa-70952_2_PmB123_2.Sam
samtools sort Pm_Solexa-70952_2_PmB123.bam Pm_Solexa-70952_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pm_Solexa-70952_2_PmB123.sorted.bam \
O=Pm_Solexa-70952_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pm_Solexa-70953
$LASTDIR-744/src/lastal -j7 -S1 -Q0 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Broad_RNAseq/Pm_Solexa-70953_Left.fasta \
| $LASTDIR-744/src/last-split \
> Pm_Solexa-70953_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q0 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Broad_RNAseq/Pm_Solexa-70953_Right.fasta \
| $LASTDIR-744/src/last-split \
> Pm_Solexa-70953_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pm_Solexa-70953_2_PmB123_1.maf > Pm_Solexa-70953_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pm_Solexa-70953_2_PmB123_2.maf > Pm_Solexa-70953_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Solexa-70953_2_PmB123_1.sam > Pm_Solexa-70953_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Solexa-70953_2_PmB123_2.sam > Pm_Solexa-70953_2_PmB123_2.Sam

samtools merge -uf  Pm_Solexa-70953_2_PmB123.bam Pm_Solexa-70953_2_PmB123_1.Sam Pm_Solexa-70953_2_PmB123_2.Sam
samtools sort Pm_Solexa-70953_2_PmB123.bam Pm_Solexa-70953_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pm_Solexa-70953_2_PmB123.sorted.bam \
O=Pm_Solexa-70953_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pm_Caspo_11
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Caspo_exp/RNA_Seq_11_2015/Caspo_2_1.fastq \
| $LASTDIR-744/src/last-split \
> Pm_Caspo_11_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Caspo_exp/RNA_Seq_11_2015/Caspo_2_2.fastq \
| $LASTDIR-744/src/last-split \
> Pm_Caspo_11_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pm_Caspo_11_2_PmB123_1.maf > Pm_Caspo_11_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pm_Caspo_11_2_PmB123_2.maf > Pm_Caspo_11_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Caspo_11_2_PmB123_1.sam > Pm_Caspo_11_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Caspo_11_2_PmB123_2.sam > Pm_Caspo_11_2_PmB123_2.Sam

samtools merge -uf  Pm_Caspo_11_2_PmB123.bam Pm_Caspo_11_2_PmB123_1.Sam Pm_Caspo_11_2_PmB123_2.Sam
samtools sort Pm_Caspo_11_2_PmB123.bam Pm_Caspo_11_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pm_Caspo_11_2_PmB123.sorted.bam \
O=Pm_Caspo_11_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pm_Caspo_8
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Caspo_exp/RNA_Seq_8_2015/Caspo_4_1.fastq \
| $LASTDIR-744/src/last-split \
> Pm_Caspo_8_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Caspo_exp/RNA_Seq_8_2015/Caspo_4_2.fastq \
| $LASTDIR-744/src/last-split \
> Pm_Caspo_8_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pm_Caspo_8_2_PmB123_1.maf > Pm_Caspo_8_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pm_Caspo_8_2_PmB123_2.maf > Pm_Caspo_8_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Caspo_8_2_PmB123_1.sam > Pm_Caspo_8_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Caspo_8_2_PmB123_2.sam > Pm_Caspo_8_2_PmB123_2.Sam

samtools merge -uf  Pm_Caspo_8_2_PmB123.bam Pm_Caspo_8_2_PmB123_1.Sam Pm_Caspo_8_2_PmB123_2.Sam
samtools sort Pm_Caspo_8_2_PmB123.bam Pm_Caspo_8_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pm_Caspo_8_2_PmB123.sorted.bam \
O=Pm_Caspo_8_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pm_control_11
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Caspo_exp/RNA_Seq_11_2015/Control_10_1.fastq \
| $LASTDIR-744/src/last-split \
> Pm_control_11_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Caspo_exp/RNA_Seq_11_2015/Control_10_2.fastq \
| $LASTDIR-744/src/last-split \
> Pm_control_11_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pm_control_11_2_PmB123_1.maf > Pm_control_11_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pm_control_11_2_PmB123_2.maf > Pm_control_11_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_control_11_2_PmB123_1.sam > Pm_control_11_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_control_11_2_PmB123_2.sam > Pm_control_11_2_PmB123_2.Sam

samtools merge -uf  Pm_control_11_2_PmB123.bam Pm_control_11_2_PmB123_1.Sam Pm_control_11_2_PmB123_2.Sam
samtools sort Pm_control_11_2_PmB123.bam Pm_control_11_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pm_control_11_2_PmB123.sorted.bam \
O=Pm_control_11_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pm_control_8
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Caspo_exp/RNA_Seq_8_2015/Control_9_1.fastq \
| $LASTDIR-744/src/last-split \
> Pm_control_8_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Caspo_exp/RNA_Seq_8_2015/Control_9_2.fastq \
| $LASTDIR-744/src/last-split \
> Pm_control_8_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pm_control_8_2_PmB123_1.maf > Pm_control_8_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pm_control_8_2_PmB123_2.maf > Pm_control_8_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_control_8_2_PmB123_1.sam > Pm_control_8_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_control_8_2_PmB123_2.sam > Pm_control_8_2_PmB123_2.Sam

samtools merge -uf  Pm_control_8_2_PmB123.bam Pm_control_8_2_PmB123_1.Sam Pm_control_8_2_PmB123_2.Sam
samtools sort Pm_control_8_2_PmB123.bam Pm_control_8_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pm_control_8_2_PmB123.sorted.bam \
O=Pm_control_8_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pm_Ca11rep3
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Caspo_exp/RNA_Seq_11_2015_rep3/Caspo_1_1.fastq \
| $LASTDIR-744/src/last-split \
> Pm_Ca11rep3_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Caspo_exp/RNA_Seq_11_2015_rep3/Caspo_1_2.fastq \
| $LASTDIR-744/src/last-split \
> Pm_Ca11rep3_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pm_Ca11rep3_2_PmB123_1.maf > Pm_Ca11rep3_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pm_Ca11rep3_2_PmB123_2.maf > Pm_Ca11rep3_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Ca11rep3_2_PmB123_1.sam > Pm_Ca11rep3_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Ca11rep3_2_PmB123_2.sam > Pm_Ca11rep3_2_PmB123_2.Sam

samtools merge -uf  Pm_Ca11rep3_2_PmB123.bam Pm_Ca11rep3_2_PmB123_1.Sam Pm_Ca11rep3_2_PmB123_2.Sam
samtools sort Pm_Ca11rep3_2_PmB123.bam Pm_Ca11rep3_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pm_Ca11rep3_2_PmB123.sorted.bam \
O=Pm_Ca11rep3_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pm_Con11rep3
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Caspo_exp/RNA_Seq_11_2015_rep3/Control_11_1.fastq \
| $LASTDIR-744/src/last-split \
> Pm_Con11rep3_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Caspo_exp/RNA_Seq_11_2015_rep3/Control_11_2.fastq \
| $LASTDIR-744/src/last-split \
> Pm_Con11rep3_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pm_Con11rep3_2_PmB123_1.maf > Pm_Con11rep3_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pm_Con11rep3_2_PmB123_2.maf > Pm_Con11rep3_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Con11rep3_2_PmB123_1.sam > Pm_Con11rep3_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Con11rep3_2_PmB123_2.sam > Pm_Con11rep3_2_PmB123_2.Sam

samtools merge -uf  Pm_Con11rep3_2_PmB123.bam Pm_Con11rep3_2_PmB123_1.Sam Pm_Con11rep3_2_PmB123_2.Sam
samtools sort Pm_Con11rep3_2_PmB123.bam Pm_Con11rep3_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pm_Con11rep3_2_PmB123.sorted.bam \
O=Pm_Con11rep3_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pm A123
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Broad_DNA-seq/A123_1.fq \
| $LASTDIR-744/src/last-split \
> Pm_A123_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Broad_DNA-seq/A123_2.fq \
| $LASTDIR-744/src/last-split \
> Pm_A123_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pm_A123_2_PmB123_1.maf > Pm_A123_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pm_A123_2_PmB123_2.maf > Pm_A123_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_A123_2_PmB123_1.sam > Pm_A123_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_A123_2_PmB123_2.sam > Pm_A123_2_PmB123_2.Sam

samtools merge -uf  Pm_A123_2_PmB123.bam Pm_A123_2_PmB123_1.Sam Pm_A123_2_PmB123_2.Sam
samtools sort Pm_A123_2_PmB123.bam Pm_A123_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pm_A123_2_PmB123.sorted.bam \
O=Pm_A123_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pm_Da1
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Broad_DNA-seq/Da1_1.fq \
| $LASTDIR-744/src/last-split \
> Pm_Da1_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Broad_DNA-seq/Da1_2.fq \
| $LASTDIR-744/src/last-split \
> Pm_Da1_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pm_Da1_2_PmB123_1.maf > Pm_Da1_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pm_Da1_2_PmB123_2.maf > Pm_Da1_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Da1_2_PmB123_1.sam > Pm_Da1_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Da1_2_PmB123_2.sam > Pm_Da1_2_PmB123_2.Sam

samtools merge -uf  Pm_Da1_2_PmB123.bam Pm_Da1_2_PmB123_1.Sam Pm_Da1_2_PmB123_2.Sam
samtools sort Pm_Da1_2_PmB123.bam Pm_Da1_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pm_Da1_2_PmB123.sorted.bam \
O=Pm_Da1_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pm_Da3
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Broad_DNA-seq/Da3_1.fq \
| $LASTDIR-744/src/last-split \
> Pm_Da3_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Broad_DNA-seq/Da3_2.fq \
| $LASTDIR-744/src/last-split \
> Pm_Da3_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pm_Da3_2_PmB123_1.maf > Pm_Da3_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pm_Da3_2_PmB123_2.maf > Pm_Da3_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Da3_2_PmB123_1.sam > Pm_Da3_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_Da3_2_PmB123_2.sam > Pm_Da3_2_PmB123_2.Sam

samtools merge -uf  Pm_Da3_2_PmB123.bam Pm_Da3_2_PmB123_1.Sam Pm_Da3_2_PmB123_2.Sam
samtools sort Pm_Da3_2_PmB123.bam Pm_Da3_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pm_Da3_2_PmB123.sorted.bam \
O=Pm_Da3_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pm_C1
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Broad_DNA-seq/C1_1.fq \
| $LASTDIR-744/src/last-split \
> Pm_C1_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Broad_DNA-seq/C1_2.fq \
| $LASTDIR-744/src/last-split \
> Pm_C1_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pm_C1_2_PmB123_1.maf > Pm_C1_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pm_C1_2_PmB123_2.maf > Pm_C1_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_C1_2_PmB123_1.sam > Pm_C1_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_C1_2_PmB123_2.sam > Pm_C1_2_PmB123_2.Sam

samtools merge -uf  Pm_C1_2_PmB123.bam Pm_C1_2_PmB123_1.Sam Pm_C1_2_PmB123_2.Sam
samtools sort Pm_C1_2_PmB123.bam Pm_C1_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pm_C1_2_PmB123.sorted.bam \
O=Pm_C1_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Pm_MS96
$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Broad_DNA-seq/MS96_1.fq \
| $LASTDIR-744/src/last-split \
> Pm_MS96_2_PmB123_1.maf

$LASTDIR-744/src/lastal -j7 -S1 -Q1 -r5 -q5 -a15 -b3 -i8 -P 4 Pmdb \
$READSDIR/Pm_Broad_DNA-seq/MS96_2.fq \
| $LASTDIR-744/src/last-split \
> Pm_MS96_2_PmB123_2.maf

$LASTDIR-744/scripts/maf-convert sam Pm_MS96_2_PmB123_1.maf > Pm_MS96_2_PmB123_1.sam
$LASTDIR-744/scripts/maf-convert sam Pm_MS96_2_PmB123_2.maf > Pm_MS96_2_PmB123_2.sam

module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_MS96_2_PmB123_1.sam > Pm_MS96_2_PmB123_1.Sam
samtools view -ht GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai Pm_MS96_2_PmB123_2.sam > Pm_MS96_2_PmB123_2.Sam

samtools merge -uf  Pm_MS96_2_PmB123.bam Pm_MS96_2_PmB123_1.Sam Pm_MS96_2_PmB123_2.Sam
samtools sort Pm_MS96_2_PmB123.bam Pm_MS96_2_PmB123.sorted

module load picard/2.1.1
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
I=Pm_MS96_2_PmB123.sorted.bam \
O=Pm_MS96_2_PmB123.sorted.dr.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

# Make the three sets of sliding windows (100 kb, 500 kb, 1 Mb) and concatenate them into a single file:
module load BEDTools/2.23.0-goolf-1.7.20
bedtools makewindows -g GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai \
-w 1000000 \
-s 100000  \
-i winnum |
awk '{print $0":1000kb"}' \
> windows_1000kb_Pm_ref.bed

bedtools makewindows -g GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai \
-w 500000 \
-s 50000  \
-i winnum |
awk '{print $0":500kb"}' \
> windows_500kb_Pm_ref.bed

bedtools makewindows -g GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta.fai \
-w 100000 \
-s 10000  \
-i winnum | \
awk '{print $0":100kb"}' \
> windows_100kb_Pm_ref.bed

cat windows_*_Pm_ref.bed > windows_Pm_ref.bed

# compute Fst w fermikit
~/utils/fermikit/htsbox/htsbox pileup -cuf \
$REFGENOME/NCBI/GCF_000349005.1_Pneumo_murina_B123_V2_genomic.fasta \
Pm_Con11rep3_2_PmB123.sorted.dr.bam \
Pm_Da3_2_PmB123.sorted.dr.bam \
Pm_B123_2_PmB123.sorted.dr.bam \
Pm_A123_2_PmB123.sorted.dr.bam \
Pm_C1_2_PmB123.sorted.dr.bam \
Pm_MS96_2_PmB123.sorted.dr.bam \
Pm_Da1_2_PmB123.sorted.dr.bam \
Pm_Caspo_8_2_PmB123.sorted.dr.bam \
Pm_control_8_2_PmB123.sorted.dr.bam \
Pm_control_11_2_PmB123.sorted.dr.bam \
Pm_Solexa-70952_2_PmB123.sorted.dr.bam \
Pm_Solexa-70951_2_PmB123.sorted.dr.bam \
Pm_Solexa-70953_2_PmB123.sorted.dr.bam \
Pm_Ca11rep3_2_PmB123.sorted.dr.bam \
Pm_Caspo_11_2_PmB123.sorted.dr.bam \
Pm_Solexa-70950_2_PmB123.sorted.dr.bam \
Pc_B80_2_PmB123.sorted.dr.bam \
Pc_Solexa-70955_2_PmB123.sorted.dr.bam \
Pc_Solexa-70956_2_PmB123.sorted.dr.bam \
Pc_Solexa-70954_2_PmB123.sorted.dr.bam \
Pc_B70_2_PmB123.sorted.dr.bam \
Pc_B50_2_PmB123.sorted.dr.bam \
Pj_SE8_2_PmB123.sorted.dr.bam \
Pj_RU7_2_PmB123.sorted.dr.bam \
Pj_RU12_2_PmB123.sorted.dr.bam \
Pj_RU817_2_PmB123.sorted.dr.bam \
Pj_RN_2_PmB123.sorted.dr.bam \
Pj_W_2_PmB123.sorted.dr.bam \
Pj_Z_2_PmB123.sorted.dr.bam \
> Pj_Pc_Pm_refPmurina.LAST_a15_b3.all_species.raw.fermikit.vcf

~/utils/fermikit/fermi.kit/k8 ~/utils/fermikit/fermi.kit/hapdip.js vcfsum \
-f Pj_Pc_Pm_refPmurina.LAST_a15_b3.all_species.raw.fermikit.vcf \
> Pj_Pc_Pm_refPmurina.LAST_a15_b3.all_species.raw.fermikit.flt.vcf 

~/utils/vcftools_0.1.13/bin/vcftools \
--vcf Pj_Pc_Pm_refPmurina.LAST_a15_b3.all_species.raw.fermikit.flt.vcf \
--recode \
--stdout | ~/utils/vcftools_0.1.13/bin/vcftools --vcf - \
--max-missing 0.5 \
--min-alleles 2 --max-alleles 6 \
--recode \
--stdout \
> Pj_Pc_Pm_refPmurina.LAST_a15_b3.all_species.fermikit.flt.optimized.vcf

~/utils/vcftools_0.1.13/bin/vcftools \
--vcf Pj_Pc_Pm_refPmurina.LAST_a15_b3.all_species.fermikit.flt.optimized.vcf \
--weir-fst-pop pop_pj_fermikit.txt --weir-fst-pop pop_pm_fermikit.txt \
--stdout | \
tail -n +2 | \
awk -F $'\t' 'BEGIN{OFS=FS}{print $1,$2-1,$2,$1":"$2,$3}' \
> pop_pj_pm_fst_refPm.fermit.optimized.bed

~/utils/vcftools_0.1.13/bin/vcftools \
--vcf Pj_Pc_Pm_refPmurina.LAST_a15_b3.all_species.fermikit.flt.optimized.vcf \
--weir-fst-pop pop_pj_fermikit.txt --weir-fst-pop pop_pc_fermikit.txt \
--stdout | \
tail -n +2 | \
awk -F $'\t' 'BEGIN{OFS=FS}{print $1,$2-1,$2,$1":"$2,$3}' \
> pop_pj_pc_fst_refPm.fermit.optimized.bed

~/utils/vcftools_0.1.13/bin/vcftools \
--vcf Pj_Pc_Pm_refPmurina.LAST_a15_b3.all_species.fermikit.flt.optimized.vcf \
--weir-fst-pop pop_pc_fermikit.txt --weir-fst-pop pop_pm_fermikit.txt \
--stdout | \
tail -n +2 | \
awk -F $'\t' 'BEGIN{OFS=FS}{print $1,$2-1,$2,$1":"$2,$3}' \
> pop_pc_pm_fst_refPm.fermit.optimized.bed

module load BEDTools/2.23.0-goolf-1.7.20
bedtools intersect \
-a windows_Pm_ref.bed \
-b pop_pj_pm_fst_refPm.fermit.optimized.bed  -wa -wb > windows_pj_pm_fst_refPm.fermit.optimized.tab

module load BEDTools/2.23.0-goolf-1.7.20
bedtools intersect \
-a windows_Pm_ref.bed \
-b pop_pj_pc_fst_refPm.fermit.optimized.bed  -wa -wb > windows_pj_pc_fst_refPm.fermit.optimized.tab

module load BEDTools/2.23.0-goolf-1.7.20
bedtools intersect \
-a windows_Pm_ref.bed \
-b pop_pc_pm_fst_refPm.fermit.optimized.bed  -wa -wb > windows_pc_pm_fst_refPm.fermit.optimized.tab

grep -v 'nan' windows_pj_pm_fst_refPm.fermit.optimized.tab > windows_pj_pm_fst_refPm.fermit.optimized.clean.tab
grep -v 'nan' windows_pj_pc_fst_refPm.fermit.optimized.tab > windows_pj_pc_fst_refPm.fermit.optimized.clean.tab
grep -v 'nan' windows_pc_pm_fst_refPm.fermit.optimized.tab > windows_pc_pm_fst_refPm.fermit.optimized.clean.tab

bedtools groupby -i windows_pj_pm_fst_refPm.fermit.optimized.clean.tab -g 1,2,3,4 -c 9 -o mean | tr "\:" "\t" > windows_pj_pm_fst_refPm.fermit.optimized.mean.tab
bedtools groupby -i windows_pj_pc_fst_refPm.fermit.optimized.clean.tab -g 1,2,3,4 -c 9 -o mean | tr "\:" "\t" > windows_pj_pc_fst_refPm.fermit.optimized.mean.tab
bedtools groupby -i windows_pc_pm_fst_refPm.fermit.optimized.clean.tab -g 1,2,3,4 -c 9 -o mean | tr "\:" "\t" > windows_pc_pm_fst_refPm.fermit.optimized.mean.tab


# PCA
~/utils/vcftools_0.1.13/bin/vcftools \
--vcf Pj_Pc_Pm_refPmurina.LAST_a15_b3.all_species.fermikit.flt.optimized.vcf \
--plink --out Pj_Pc_Pm_refPmurina.LAST_a15_b3.all_species.fermikit.flt.optimized.vcf.raw

#Make a .genome file
~/utils/plink --file Pj_Pc_Pm_refPmurina.LAST_a15_b3.all_species.fermikit.flt.optimized.vcf.raw --genome --noweb --allow-no-sex \
--out Pj_Pc_Pm_refPmurina.LAST_a15_b3.all_species.fermikit.flt.optimized.vcf.raw

#Do some multidimensional scaling:
~/utils/plink --file Pj_Pc_Pm_refPmurina.LAST_a15_b3.all_species.fermikit.flt.optimized.vcf.raw \
-read-genome Pj_Pc_Pm_refPmurina.LAST_a15_b3.all_species.fermikit.flt.optimized.vcf.raw.genome --cluster --mds-plot 2 --noweb


# plotting in R
fst <- read.table("$WRKDIR/windows_pj_pc_fst_refPm.fermit.optimized.mean.tab",header=F,sep="\t")
names(fst) <- c("chrom", "start", "end", "win_id","win_size", "avg_fst" )

fst$win_size <- factor(fst$win_size, levels=c("100kb", "500kb", "1000kb"))

ggplot(data=fst, aes(x=avg_fst)) +
  geom_density(fill=I("blue")) +
  facet_wrap(~win_size)

ggplot(fst, aes(y=avg_fst, x=start, colour=win_size)) +
  geom_line() +
  facet_wrap(~chrom, nrow=2) +
  scale_colour_manual(name="Window size", values=c("green", "blue","red"))

q <- quantile(subset(fst,win_size=="500kb",select="avg_fst")[,1],prob=0.99)[[1]]

ggplot(fst, aes(y=avg_fst, x=start, colour=win_size)) +
  geom_line() +
  facet_wrap(~chrom, nrow=2) +
  geom_hline(yintercept=q,colour="black") +
  scale_colour_manual(name="Window size", values=c("green", "blue","red"))

previous_theme <- theme_set(theme_bw())

ggplot(data=fst, aes(x=avg_fst)) +
  geom_density(fill=I("blue")) +
  facet_wrap(~win_size) +  geom_vline(xintercept=q,colour="black")

d<-read.table("plink.txt",h=T)
plot(d$C1,d$C2,  col=as.integer(d$FID), pch=19, xlab = "", ylab = "", axes = T, main = "PCA of Pneumocystis samples/plink (SNPs)")
legend("topright", c("P.jirovecii", "P.carinii","P.murina"), pch=19, col=c(2,1,3))
text(d$C1,d$C2, labels =  (d$IID), pos = 4)

