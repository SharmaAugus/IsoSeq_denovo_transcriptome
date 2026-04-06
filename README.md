# MORE24 Flower Iso-Seq Transcriptome Project

## Overview

This repository documents the construction of a **de novo reference transcriptome** for a flower species (*Brassicales*) using Pacific Biosciences **Iso-Seq long-read sequencing** data, in the absence of a reference genome.

---

## Key Results

| Metric | Value |
|--------|-------|
| Final reference sequences | 82,994 |
| N50 | 1,979 bp |
| BUSCO completeness (Brassicales) | 95.3% |
| Complete ORFs (TransDecoder) | 85.2% |
| Gene families identified | 27,219 |

---

## Pipeline Summary

```
Raw FLNC Iso-Seq reads (2 replicates)
        ↓
CD-HIT-EST clustering (95% identity, >500bp)
        ↓
Cogent: K-mer profiling + Gene family finding
        ↓
Cogent: Coding genome reconstruction
        ↓
Rename + CD-HIT-EST final cleanup
        ↓
MORE24_FL_final_reference_nr.fa (82,994 sequences)
        ↓
BUSCO + TransDecoder quality assessment
```

---

## Repository Contents

| File/Folder | Description |
|-------------|-------------|
| `README.md` | This file — project overview |
| `docs/De_novo_Reference_Transcriptome_Cogent_Pipeline.sh` | Full pipeline with all commands |
| `scripts/rename_final_reference_v2.py` | Sequence renaming script |
| `scripts/Length_Dist_final.Ref.py` | Sequence length statistics script |

---

## Tools Used

[Cogent](https://github.com/Magdoll/Cogent) •
[cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake) •
[CD-HIT](https://github.com/weizhongli/cdhit) •
[minimap2](https://github.com/lh3/minimap2) •
[MASH](https://github.com/marbl/Mash) •
[BUSCO](https://busco.ezlab.org/) •
[TransDecoder](https://github.com/TransDecoder/TransDecoder)

---

## Documentation

For the complete step-by-step pipeline with all commands, parameters, and results, see:

📄 **[De novo Reference Transcriptome Pipeline Using Cogent](docs/De_novo_Reference_Transcriptome_Cogent_Pipeline.md)**

---

## Contact at saloni@ugr.go.es for any information

**Project:** MORE24 — Flower Iso-Seq Transcriptome
**Species:** Brassicales (flower, no reference genome)
**Data:** Pacific Biosciences Iso-Seq (2 biological replicates)
