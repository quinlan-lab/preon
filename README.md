# preon
Precision Oncology data workflow

CADD downloads:
+ wget -c https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz
+ wget -c https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz.tbi
+ wget -c https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/InDels.tsv.gz
+ wget -c https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/InDels.tsv.gz.tbi

Manifest file:

```
$ cat manifest.tsv
tow19   tow19_normal    0       190719_A00421_0090_BHKNYWDMXX_S1        /uufs/chpc.utah.edu/common/HIPAA/u6022494/work/bwa-nf/test-data/16239X1_190719_A00421_0090_BHKNYWDMXX_S1_L002_R1_001.fastq.gz    /uufs/chpc.utah.edu/common/HIPAA/u6022494/work/bwa-nf/test-data/16239X1_190719_A00421_0090_BHKNYWDMXX_S1_L002_R2_001.fastq.gz
tow19   tow19_tumor     1       190628_A00421_0082_AHKTTKDMXX_S19       /uufs/chpc.utah.edu/common/HIPAA/u6022494/work/bwa-nf/test-data/16201X1_190628_A00421_0082_AHKTTKDMXX_S19_L001_R1_001.fastq.gz   /uufs/chpc.utah.edu/common/HIPAA/u6022494/work/bwa-nf/test-data/16201X1_190628_A00421_0082_AHKTTKDMXX_S19_L001_R2_001.fastq.gz
tow19   tow19_tumor     1       190628_A00421_0082_AHKTTKDMXX_S20       /uufs/chpc.utah.edu/common/HIPAA/u6022494/work/bwa-nf/test-data/16201X2_190628_A00421_0082_AHKTTKDMXX_S20_L001_R1_001.fastq.gz   /uufs/chpc.utah.edu/common/HIPAA/u6022494/work/bwa-nf/test-data/16201X2_190628_A00421_0082_AHKTTKDMXX_S20_L001_R2_001.fastq.gz
```

Usage:

```
nextflow run main.nf -c nextflow.config -resume -profile redwood -qs 100 \
    --project tow19 --manifest manifest.tsv \
    --fasta /uufs/chpc.utah.edu/common/HIPAA/u6022494/work/assets/broad_bundle_hg38/Homo_sapiens_assembly38.fasta \
    --caddsnv /uufs/chpc.utah.edu/common/HIPAA/u6022494/work/assets/GRCh38/cadd/whole_genome_SNVs.tsv.gz \
    --caddindel /uufs/chpc.utah.edu/common/HIPAA/u6022494/work/assets/GRCh38/cadd/InDels.tsv.gz
```
