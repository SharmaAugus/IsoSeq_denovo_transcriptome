from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

output_records = []
gene_counter = 1
total_written = 0

print("Processing reconstructed families...")

# Step 1: Reconstructed families - read directly, don't use SeqIO on combined file
family_dirs = sorted(os.listdir('flower_output'))

for family_name in family_dirs:
    family_path = os.path.join('flower_output', family_name)
    cogent_fa = os.path.join(family_path, 'cogent2.fa')
    
    if not os.path.exists(cogent_fa):
        continue
    
    isoform_counter = 0
    for rec in SeqIO.parse(cogent_fa, 'fasta'):
        new_id = f"MORE24_FL_gene{gene_counter}_isoform{gene_counter}.{isoform_counter}"
        new_rec = SeqRecord(
            rec.seq,
            id=new_id,
            description=f"original_family={family_name}"
        )
        output_records.append(new_rec)
        isoform_counter += 1
        total_written += 1
    
    gene_counter += 1

print(f"Reconstructed: {total_written} sequences from {gene_counter-1} families")
recon_count = total_written

# Step 2: Unresolved families
print("\nProcessing unresolved families...")
unresolved_families = ['flower_106', 'flower_128', 'flower_25', 'flower_6']

unresolved_count = 0
for family_name in unresolved_families:
    in_fa = os.path.join('flower_output', family_name, 'in.fa')
    if not os.path.exists(in_fa):
        continue
    
    isoform_counter = 0
    for rec in SeqIO.parse(in_fa, 'fasta'):
        new_id = f"MORE24_FL_gene{gene_counter}_isoform{gene_counter}.{isoform_counter}"
        new_rec = SeqRecord(
            rec.seq,
            id=new_id,
            description=f"original_family={family_name} unresolved original_id={rec.id}"
        )
        output_records.append(new_rec)
        isoform_counter += 1
        total_written += 1
        unresolved_count += 1
    
    gene_counter += 1

print(f"Unresolved: {unresolved_count} sequences")

# Step 3: Singletons
print("\nProcessing singletons...")

singleton_ids = set()
with open('flower.partition.txt') as f:
    for line in f:
        if line.startswith('#unassigned:'):
            members = line.strip().replace('#unassigned:', '')
            for member in members.split(','):
                singleton_ids.add(member.strip())

print(f"Singleton IDs loaded: {len(singleton_ids)}")

singleton_count = 0
for rec in SeqIO.parse('Fl.isoseq_flnc.fasta', 'fasta'):
    if rec.id in singleton_ids:
        new_id = f"MORE24_FL_gene{gene_counter}_isoform{gene_counter}.0"
        new_rec = SeqRecord(
            rec.seq,
            id=new_id,
            description=f"original_id={rec.id} singleton"
        )
        output_records.append(new_rec)
        gene_counter += 1
        total_written += 1
        singleton_count += 1

print(f"Singletons: {singleton_count} sequences")

# Write output - write one by one to avoid any deduplication
print(f"\nWriting output...")
print(f"Total genes:     {gene_counter-1}")
print(f"Total sequences: {total_written}")
print(f"  Reconstructed: {recon_count}")
print(f"  Unresolved:    {unresolved_count}")
print(f"  Singletons:    {singleton_count}")

with open('MORE24_FL_final_reference_v2.fa', 'w') as out:
    for rec in output_records:
        out.write(f">{rec.id} {rec.description}\n")
        # Write sequence in 60 char lines
        seq = str(rec.seq)
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + '\n')

# Verify
count = sum(1 for line in open('MORE24_FL_final_reference_v2.fa') if line.startswith('>'))
print(f"\nVerification count: {count} sequences")

# Sample headers
print("\nSample headers:")
with open('MORE24_FL_final_reference_v2.fa') as f:
    shown = 0
    for line in f:
        if line.startswith('>'):
            print(f"  {line.strip()}")
            shown += 1
            if shown >= 5:
                break
