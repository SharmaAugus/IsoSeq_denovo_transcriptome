# De novo Reference Transcriptome Pipeline Using Cogent
### MORE24 Flower Iso-Seq Project

## Project Overview

Previous 'de novo reference transcriptome' for a flower species (Moricandia. arvensis) without a reference genome, using Pacific Biosciences Iso-Seq long-read data, reclustered with CD-HIT has some high number of isoforms per gene (1000s).
Thus here, using  'Cogent' coding genome reconstruction final de novo reference transcriptome is generated, and quality assessment performed as previous pipeline by BUSCO & TransDecoder.
---

## Final Results Summary

| Metric | Value |
|--------|-------|
| Final reference sequences | 82,994 |
| Min sequence length | 500 bp |
| Max sequence length | 17,572 bp |
| Mean length | 1,828 bp |
| N50 | 1,979 bp |
| Total bases | ~152 Mb |
| BUSCO completeness (Brassicales) | 95.3% |
| Complete ORFs | 85.2% |

---

## Pipeline Overview

```
Raw FLNC reads (1.15M)
    ↓ Iso-Seq Clustering
Clustered HQ transcripts (0.8M)
    ↓ CD-HIT-EST (95% identity, >500bp)
Non-redundant input (160,410)
    ↓ Cogent: run_mash.py
    ↓ Cogent: process_kmer_to_graph.py
27,219 gene families
    ↓ Cogent: reconstruct_contig.py
86,618 reconstructed contigs
    + 126 unresolved sequences
    + 6,857 singletons
= 93,601 combined
    ↓ Rename headers
    ↓ CD-HIT-EST (95% final cleanup)
82,994 Final Reference Transcriptome
    ↓ BUSCO (Brassicales)       → 95.3% complete
    ↓ TransDecoder              → 85.2% complete ORFs
```

---

## Software Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| Python | 3.7 | Core scripting |
| Conda | Any | Environment management |
| cDNA_Cupcake | 29.0.0 | Iso-Seq post-processing |
| Cogent | 8.0.0 | Coding genome reconstruction |
| CD-HIT-EST | Any | Sequence clustering |
| minimap2 | 2.28 | Sequence alignment |
| mash | 1.1 | K-mer distance calculation |
| BUSCO | 5.x | Transcriptome completeness |
| TransDecoder | Any | ORF prediction |
| BioPython | Any | Sequence parsing |
| parasail | Any | Sequence alignment (Cogent dep) |
| PuLP | Any | LP solver (Cogent dep) |

---

## Step 1: Environment Setup

```bash
# Create conda environment
conda create -n cogent_env python=3.7 -y
conda activate cogent_env

# Install Python dependencies
pip install cython
pip install biopython bcbio-gff networkx pulp parasail scikit-learn

# Install mash
conda install -c bioconda mash -y

# Install minimap2
conda install -c bioconda minimap2 -y

# Install GNU parallel
conda install -c conda-forge parallel -y
```

---

## Step 2: Install cDNA_Cupcake

```bash
git clone https://github.com/Magdoll/cDNA_Cupcake.git
cd cDNA_Cupcake
pip install -e .
cd ..
```

---

## Step 3: Install Cogent

```bash
git clone https://github.com/Magdoll/Cogent.git
cd Cogent
pip install -e .
cd ..

# Set path variable for convenience
export COGENT=/path/to/Cogent/Cogent
```

---

## Step 4: Input Preparation

### 4.1 Combine FLNC reads from replicates

Two biological replicates were combined before clustering:

```bash
cat Replicate1_FLNC.fasta Replicate2_FLNC.fasta > combined_FLNC.fasta
```
# Here used isoseq clustered hq.fa file (file from previous command)


### 4.2 CD-HIT-EST clustering (95% identity, >500bp filter)

```bash
cd-hit-est \
    -i combined_FLNC.fasta \
    -o Flower.nonredundant.fasta \
    -c 0.95 -n 8 \
    -T 20 -M 0 \
    -d 0 \
    -min_l 500 \
    > cdhit_input.log 2>&1

# Result: 160,410 sequences
grep -c ">" Flower.nonredundant.fasta
```

---

## Step 5: Cogent Analysis

### 5.1 Setup working directory

```bash
mkdir Cogent_Analysis
cd Cogent_Analysis

# Symlink input file (required by Cogent)
ln -s /path/to/Flower.nonredundant.fasta Fl.isoseq_flnc.fasta

# Verify
grep -c ">" Fl.isoseq_flnc.fasta   # should show 160,410
```

### 5.2 K-mer profiling with MASH

```bash
# This step generates pairwise k-mer distances
# Expected output: Fl.isoseq_flnc.fasta.s1000k30.dist (~48G)
# Expected runtime: 3-6 hours

nohup run_mash.py -k 30 --cpus=20 Fl.isoseq_flnc.fasta \
    > run_mash.log 2>&1 &

# Monitor progress
grep "completed" run_mash.log | wc -l
```

### 5.3 Graph-based family finding

```bash
# Partitions 160,410 sequences into gene families
# Expected output: flower_output/ directory + flower.partition.txt
# Expected runtime: 2-4 hours

nohup process_kmer_to_graph.py \
    Fl.isoseq_flnc.fasta \
    Fl.isoseq_flnc.fasta.s1000k30.dist \
    flower_output/ \
    flower \
    > process_kmer.log 2>&1 &

# Check results when complete
wc -l flower.partition.txt          # number of gene families
ls flower_output/ | wc -l           # number of output directories
```

**Result:** 27,219 gene families identified

```bash
# Check family size distribution
awk -F'\t' 'NR>1{print $2}' flower.partition.txt | sort -n | uniq -c

# Count unassigned sequences
grep "unassigned" flower.partition.txt | grep -o "transcript/[0-9]*" | wc -l
# Result: 6,857 singletons

# Largest gene families
awk -F'\t' 'NR>1{print $2, $1}' flower.partition.txt | sort -rn | head -10
```

### 5.4 Coding genome reconstruction

```bash
# Generate batch commands for all families
nohup generate_batch_cmd_for_Cogent_reconstruction.py \
    flower_output/ \
    > generate_reconstruction.log 2>&1 &

# Extract commands to file
grep "reconstruct_contig.py" generate_reconstruction.log \
    > reconstruction_cmds.txt

wc -l reconstruction_cmds.txt   # should show 27,218

# Run all reconstruction jobs in parallel (100 simultaneous)
# Expected runtime: 3-8 hours depending on server
nohup parallel -j 100 < reconstruction_cmds.txt \
    > reconstruction.log 2>&1 &

# Monitor progress
watch -n 60 "echo 'Completed:' && \
    ls flower_output/*/cogent2.fa 2>/dev/null | wc -l && \
    echo 'Total: 27218'"
```

### 5.5 Handle failed reconstructions

```bash
# Find failed families
for d in flower_output/flower_*/; do
    if [ ! -f "$d/cogent2.fa" ]; then
        echo "$d"
    fi
done

# Retry with cycle detection and larger k-mer
# (replace family names with your actual failed families)
cat > rerun_failed.txt << 'EOF'
reconstruct_contig.py --nx_cycle_detection -k 40 flower_output/flower_XXX -p flower_XXX
EOF

nohup parallel -j 18 < rerun_failed.txt > rerun_failed.log 2>&1 &

# Final check - count successes
ls flower_output/*/cogent2.fa 2>/dev/null | wc -l
# Result: 27,214 successful (4 permanently failed)
```

---

## Step 6: Build Final Reference Transcriptome

### 6.1 Combine all reconstructed sequences

```bash
# Combine cogent2.fa from all successful families
cat flower_output/*/cogent2.fa > cogent_reconstructed.fa
grep -c ">" cogent_reconstructed.fa   # 86,618

# Extract raw sequences from 4 permanently failed families
for d in flower_output/flower_106 flower_output/flower_128 \
          flower_output/flower_25 flower_output/flower_6; do
    cat $d/in.fa >> unresolved_families_raw.fa
done

# Deduplicate unresolved sequences
python3 -c "
from Bio import SeqIO
seen = set()
records = []
for rec in SeqIO.parse('unresolved_families_raw.fa', 'fasta'):
    if rec.id not in seen:
        seen.add(rec.id)
        records.append(rec)
with open('unresolved_families.fa', 'w') as out:
    SeqIO.write(records, out, 'fasta')
print(f'Unique unresolved sequences: {len(records)}')
"
# Result: 126 unique sequences
```

### 6.2 Extract singleton sequences

```bash
python3 -c "
from Bio import SeqIO

# Load singleton IDs from partition file
ids = set()
with open('flower.partition.txt') as f:
    for line in f:
        if line.startswith('#unassigned:'):
            members = line.strip().replace('#unassigned:', '')
            for member in members.split(','):
                ids.add(member.strip())

print(f'Singleton IDs loaded: {len(ids)}')

# Extract from original fasta
count = 0
with open('singleton_sequences.fa', 'w') as out:
    for rec in SeqIO.parse('Fl.isoseq_flnc.fasta', 'fasta'):
        if rec.id in ids:
            SeqIO.write(rec, out, 'fasta')
            count += 1
print(f'Extracted: {count} sequences')
"
# Result: 6,857 singleton sequences
```

### 6.3 Rename sequences with meaningful IDs

```bash
# Creates MORE24_FL_gene{N}_isoform{N}.{M} naming scheme
python3 rename_final_reference_v2.py

# Clean duplicate description from headers
sed 's/\(>MORE24_FL_gene[^ ]*\) MORE24_FL_gene[^ ]* /\1 /' \
    MORE24_FL_final_reference_v2.fa \
    > MORE24_FL_final_reference_v2_clean.fa

# Remove spaces in description for compatibility
sed 's/ original_family=.*$//' \
    MORE24_FL_final_reference_v2_clean.fa \
    > MORE24_FL_final_reference_clean.fa

# Verify
grep -c ">" MORE24_FL_final_reference_clean.fa   # 93,601
grep ">" MORE24_FL_final_reference_clean.fa | head -5
```

### 6.4 Final CD-HIT-EST clustering

```bash
# Remove any new redundancy introduced during reconstruction
nohup cd-hit-est \
    -i MORE24_FL_final_reference_clean.fa \
    -o MORE24_FL_final_reference_nr.fa \
    -c 0.95 -n 8 \
    -T 32 -M 0 -d 0 \
    > cdhit_final.log 2>&1 &

# Final count
grep -c ">" MORE24_FL_final_reference_nr.fa   # 82,994
```

---

## Step 7: Quality Assessment

### 7.1 Sequence statistics

```python
# save as: Length_Dist_final.Ref.py
from Bio import SeqIO

lengths = sorted([len(r) for r in SeqIO.parse('MORE24_FL_final_reference_nr.fa', 'fasta')])
total = len(lengths)
total_bases = sum(lengths)

print(f'Total sequences:  {total}')
print(f'Min length:       {min(lengths)} bp')
print(f'Max length:       {max(lengths)} bp')
print(f'Mean length:      {sum(lengths)/total:.0f} bp')
print(f'Median length:    {lengths[total//2]} bp')
print(f'Total bases:      {total_bases:,} bp')

# N50 calculation
cumsum = 0
for l in sorted(lengths, reverse=True):
    cumsum += l
    if cumsum >= total_bases / 2:
        print(f'N50:              {l} bp')
        break
```

```bash
python Length_Dist_final.Ref.py
```

**Results:**
```
Total sequences:  82,994
Min length:       500 bp
Max length:       17,572 bp
Mean length:      1,828 bp
Median length:    1,673 bp
N50:              1,979 bp
Total bases:      151,724,884 bp
```

### 7.2 BUSCO assessment

```bash
conda activate busco5

busco \
    -i MORE24_FL_final_reference_nr.fa \
    -l /path/to/busco_downloads/lineages/brassicales_odb10 \
    -o busco_brassicales \
    -m transcriptome \
    -c 20
```

**Results:**
```
C:95.3%[S:37.8%,D:57.5%],F:0.9%,M:3.9%,n:4596
Complete BUSCOs (C):            4,378  (95.3%)
  Single-copy (S):              1,735  (37.8%)
  Duplicated (D):               2,643  (57.5%)
Fragmented BUSCOs (F):             40   (0.9%)
Missing BUSCOs (M):               178   (3.9%)
Total BUSCO groups searched:    4,596
```

> Note: High duplication (57.5%) is expected for transcriptomes due to multiple isoforms per gene.

---

## Step 8: Isoform Analysis

### 8.1 Isoforms per gene (Cogent-based)

```bash
python3 -c "
import os

with open('cogent_isoforms_per_gene.txt', 'w') as out:
    out.write('Gene_Family\tInput_Transcripts\tOutput_Isoforms\n')
    for family in sorted(os.listdir('flower_output')):
        family_path = os.path.join('flower_output', family)
        in_fa = os.path.join(family_path, 'in.fa')
        cogent_fa = os.path.join(family_path, 'cogent2.fa')
        if not os.path.exists(in_fa):
            continue
        input_count = sum(1 for l in open(in_fa) if l.startswith('>'))
        output_count = sum(1 for l in open(cogent_fa) if l.startswith('>')) \
                      if os.path.exists(cogent_fa) else 0
        out.write(f'{family}\t{input_count}\t{output_count}\n')
print('Done: cogent_isoforms_per_gene.txt')
"
```

### 8.2 Isoform summary table

```bash
awk '
NR>1 {
    total++
    n=$3
    if (n == 0) zero++
    else if (n == 1) single++
    else if (n >= 2 && n <= 5) typical++
    else if (n >= 6 && n <= 10) moderate++
    else if (n >= 11 && n <= 20) complex_low++
    else if (n >= 21 && n <= 50) complex_mid++
    else if (n >= 51 && n <= 100) complex_high++
    else if (n > 100) suspicious++
}
END {
    print "Category\tIsoform_Range\tCount\tPercent"
    printf "Unresolved\t0\t%d\t%.2f%%\n", zero, zero/total*100
    printf "Single\t1\t%d\t%.2f%%\n", single, single/total*100
    printf "Typical\t2-5\t%d\t%.2f%%\n", typical, typical/total*100
    printf "Moderate\t6-10\t%d\t%.2f%%\n", moderate, moderate/total*100
    printf "Complex_low\t11-20\t%d\t%.2f%%\n", complex_low, complex_low/total*100
    printf "Complex_mid\t21-50\t%d\t%.2f%%\n", complex_mid, complex_mid/total*100
    printf "Complex_high\t51-100\t%d\t%.2f%%\n", complex_high, complex_high/total*100
    printf "Suspicious\t>100\t%d\t%.2f%%\n", suspicious, suspicious/total*100
}
' cogent_isoforms_per_gene.txt > FL_cogent_isoform_summary.tsv

cat FL_cogent_isoform_summary.tsv
```

**Results:**
```
Category      Isoform_Range  Count   Percent
Unresolved    0              4       0.01%
Single        1              6,141   22.56%
Typical       2-5            17,575  64.57%
Moderate      6-10           2,963   10.89%
Complex_low   11-20          487     1.79%
Complex_mid   21-50          44      0.16%
Complex_high  51-100         3       0.01%
Suspicious    >100           1       0.00%
```

---

## Step 9: ORF Prediction with TransDecoder

```bash
conda activate cogent_env

# Step 1: Extract long ORFs
TransDecoder.LongOrfs -t MORE24_FL_final_reference_nr.fa

# Step 2: Predict coding regions
TransDecoder.Predict \
    -t MORE24_FL_final_reference_nr.fa \
    --single_best_only

# Step 3: Check ORF type distribution
grep ">" MORE24_FL_final_reference_nr.fa.transdecoder.pep | \
    grep -o "type:[^ ]*" | sort | uniq -c | sort -rn
```

**Results:**
```
60,585  type:complete         (85.2%)
 9,276  type:5prime_partial   (13.0%)
 1,197  type:3prime_partial    (1.7%)
    65  type:internal          (0.1%)
─────────────────────────────────────
71,123  Total predicted ORFs
```

---

## Output Files

| File | Description |
|------|-------------|
| `MORE24_FL_final_reference_nr.fa` | **Final reference transcriptome** (82,994 sequences) |
| `MORE24_FL_final_reference_nr.fa.transdecoder.pep` | Predicted protein sequences |
| `MORE24_FL_final_reference_nr.fa.transdecoder.cds` | Predicted CDS sequences |
| `MORE24_FL_final_reference_nr.fa.transdecoder.gff3` | ORF coordinates |
| `flower.partition.txt` | Gene family assignments |
| `cogent_isoforms_per_gene.txt` | Isoform counts per gene family |
| `FL_cogent_isoform_summary.tsv` | Isoform distribution summary |
| `singleton_sequences.fa` | Unassigned singleton sequences (6,857) |
| `busco_brassicales/` | BUSCO assessment results |

---

## Key Helper Scripts

### rename_final_reference_v2.py
Renames sequences from Cogent path-based IDs to meaningful `MORE24_FL_gene{N}_isoform{N}.{M}` format.

### Length_Dist_final.Ref.py
Calculates sequence length statistics including N50.

---

## Sequence Accounting

```
Input to Cogent:              160,410
  Assigned to families:   153,553
  Singletons:               6,857

Gene families:                 27,219
  Successfully reconstructed: 27,214
  Failed (unresolvable cycles):    4 (flower_106, flower_128, flower_25, flower_6)

Final reference components:
   Reconstructed contigs:   86,618
   Unresolved raw seqs:        126
  Singletons:               6,857
  = Combined:                  93,601

 After CD-HIT-EST 95%:      82,994 (final)
```

---

## Citation

If you use this pipeline, please cite:

- **Cogent**: Magdoll, E. (2020). Cogent v8.0.0. https://github.com/Magdoll/Cogent
- **cDNA_Cupcake**: Magdoll, E. https://github.com/Magdoll/cDNA_Cupcake
- **CD-HIT**: Fu L. et al. (2012) CD-HIT: accelerated for clustering the next-generation sequencing data. Bioinformatics.
- **MASH**: Ondov BD et al. (2016) Mash: fast genome and metagenome distance estimation using MinHash. Genome Biology.
- **minimap2**: Li H. (2018) Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics.
- **BUSCO**: Manni M. et al. (2021) BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage. Molecular Biology and Evolution.
- **TransDecoder**: Haas B. https://github.com/TransDecoder/TransDecoder

---

## Contact saloni@ugr.go.es

Project: MORE24 — Flower Iso-Seq Transcriptome
