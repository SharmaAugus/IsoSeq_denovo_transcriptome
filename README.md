# IsoSeq Full Pipeline

## Overview
Complete IsoSeq analysis workflow including:
- FLNC processing
- Clustering
- Annotation
- ORF prediction
- Quantification

## Requirements
- samtools
- isoseq3
- minimap2
- seqkit
- blast+
- cd-hit
- TransDecoder
- busco

## Run pipeline
```bash
chmod +x isoseq_pipeline.sh
./isoseq_pipeline.sh
