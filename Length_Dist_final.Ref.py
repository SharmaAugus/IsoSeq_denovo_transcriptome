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
