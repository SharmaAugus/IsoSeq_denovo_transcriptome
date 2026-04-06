#!/bin/bash
# https://isoseq.how/clustering/cli-workflow.html

############################################
# CONFIGURATION
############################################
THREADS=8
PRIMERS="Primers.fasta"
BLAST_DB="/home/saloni/blastdb/nr"

#Step 1 - Circular Consensus Sequence calling

############################################
# STEP 2: LIMA (Primer removal)
############################################
# If there are two sequences in your `Primer.fasta` file

# Primer fasta of Iso-seq:
# >Isoseq_5p
# GCAATGAAGTCGCAGGGTTGGG
# >IsoSeq_3p
# GTACTCTGCGTTGATACCACTGCTT

############################################
echo "Running LIMA..."

lima MORE24-FP.HiFi_reads.bam Primers.fasta  Isoseq_output/MORE24.FP.flnc.bam --isoseq
lima MORE24-FV.HiFi_reads.bam Primers.fasta  Isoseq_output/MORE24.FV.flnc.bam --isoseq
lima MORE24-HP.HiFi_reads.bam Primers.fasta  Isoseq_output/MORE24.HP.flnc.bam --isoseq
lima MORE24-HV.HiFi_reads.bam Primers.fasta  Isoseq_output/MORE24.HV.flnc.bam --isoseq

############################################
# STEP 3: REFINE
############################################
echo "Running IsoSeq refine..."

isoseq refine MORE24.FP.flnc.Isoseq_5p--IsoSeq_3p.bam Primers.fasta Isoseq_output/MORE24-FP.flnc.bam
isoseq refine MORE24.FV.flnc.Isoseq_5p--IsoSeq_3p.bam Primers.fasta Isoseq_output/MORE24-FV.flnc.bam
isoseq refine MORE24.HP.flnc.Isoseq_5p--IsoSeq_3p.bam Primers.fasta Isoseq_output/MORE24-HP.flnc.bam
isoseq refine MORE24.HV.flnc.Isoseq_5p--IsoSeq_3p.bam Primers.fasta Isoseq_output/MORE24-HV.flnc.bam

############################################
# STEP 4: MERGE FLNC FILES
############################################
echo "Merging FLNC BAM files..."

samtools merge All_flnc.bam \
  MORE24-FP.flnc.bam \
  MORE24.FV.flnc.bam \
  MORE24-HP.flnc.bam \
  MORE24-HV.flnc.bam


############################################
# STEP 5: CLEAN BAM (REMOVE RG TAGS)
############################################
echo "Cleaning BAM (removing RG tags)..."

samtools view -h All_flnc.bam > All_flnc.full.sam
grep -v "^@RG" All_flnc.full.sam | sed 's/\tRG:Z:[^\t]*//g' > All_flnc.noRG.sam
samtools view -bS All_flnc.noRG.sam > All_flnc.noRG.bam

# same for merging Flower and leaf samples separately
samtools view -h Flower_flnc.bam > Flower_flnc.full.sam
grep -v "^@RG" Flower_flnc.full.sam | sed 's/\tRG:Z:[^\t]*//g' > Flower_flnc.noRG.sam
samtools view -bS Flower_flnc.noRG.sam > Flower_flnc.noRG.bam
isoseq3 cluster Flower_flnc.noRG.bam Cluster/Flower.clustered.bam

############################################
# STEP 6: CLUSTER
############################################
echo "Running IsoSeq clustering..."

isoseq3 cluster All_flnc.noRG.bam Cluster/All.clustered.bam

# Similarly for Flower and leaf samples.

############################################
# STEP 7: EXTRACT TRANSCRIPTS
############################################
echo "Extracting transcript subsets..."

seqkit grep -f C_transcript.ids.txt All_clustered.fasta > C_extracted_transcripts.fasta
seqkit grep -f F_transcript.ids.txt MORE24.FP.clustered.hq.fasta > F_extracted_transcripts.fasta
seqkit grep -f H_transcript.ids.txt MORE24.HP.clustered.hq.fasta > H_extracted_transcripts.fasta

############################################
# STEP 8: BLASTX ANNOTATION
############################################
echo "Running BLASTX..."

/usr/bin/blastx \
  -query C_extracted_transcripts.fasta \
  -db $BLAST_DB \
  -out Annotation/C_annotation-results.txt \
  -evalue 1e-5 \
  -outfmt 6 \
  -max_target_seqs 10 \
  -num_threads 4


############################################
# STEP 9: MINIMAP2 (SELF-MAPPING)
############################################
echo "Running minimap2..."

minimap2 -ax splice -t $THREADS \
  Cluster/All.clustered.hq.fasta \
  Cluster/All.clustered.hq.fasta > mapped.sam

############################################
# STEP 10: TAMA COLLAPSE
############################################
echo "Running TAMA collapse..."

python tama_collapse.py \
  -s mapped.sam \
  -f Cluster/All.clustered.hq.fasta \
  -p TAMA/output

############################################
# STEP 11: FILTERING
############################################
echo "Filtering transcripts..."

awk '$7 != 0' TAMA/output.bed > step1.bed
awk '$7 >= 100' step1.bed > step2.bed
awk '$10 >= 3' step2.bed > step3.bed
awk '($3 - $2) >= 300' step3.bed > final_filtered.bed



# For De novo Reference use Step 6 outputs
############################################

# STEP 12: CD-HIT (REMOVE REDUNDANCY)
############################################
echo "Running CD-HIT..."

awk '/^>/ {print $0; getline seq; if(length(seq) >= 500) print seq}' /scratch/MORE24/Secuencias_Isoseq/Isoseq_output/Cluster/Flower.clustered.hq.fasta > Flower.filtered.fasta

cd-hit-est \
  -i Flower.filtered.fasta \
  -o Flower.nonredundant.fasta \
  -c 0.95 -n 10 -T 4 -M 8000

# Same for leaf samples

############################################
# STEP 13: TRANSDECODER (ORF PREDICTION)
############################################
echo "Running TransDecoder..."

TransDecoder.LongOrfs -t Flower.nonredundant.fasta
TransDecoder.Predict -t Flower.nonredundant.fasta

# Same for leaf samples


# STEP 14: BUSCO
############################################
echo "Running BUSCO..."

# cleaning of headers
awk '/^>/ {print ">transcript_" ++i} !/^>/' Flower.nonredundant.fasta > Flower.nonredundant_cleaned.cds

busco -i Flower.nonredundant_cleaned.fasta \
  -l brassicales_odb10 \
  -o busco_output \
  -m transcriptome \
  -f


# Same for leaf samples

############################################

# RSEM and Salmon  used for reference generated in both cases  through R
############################################