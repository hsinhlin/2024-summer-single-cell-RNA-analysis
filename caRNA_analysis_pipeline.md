# Chromatin-Associated RNA Sequencing (caRNA-seq) Analysis Pipeline

## Introduction
This repository documents the analysis pipeline used to process chromatin-associated RNA (caRNA-seq) data from human A375 cells (shControl vs. shMETTL3 knockdown).


The pipeline is based on the experimental and computational procedures described in: 

1. Kang, Z., Li, R., Liu, C., … Yang, X., Liu, J. **m6A‑modified cenRNA stabilizes CENPA to ensure centromere integrity in cancer cells.** *Cell* (2024). https://www.cell.com/cell/fulltext/S0092-8674(24)00969-3
2. Liu, J., Dou, X., Chen, C., Chen, C., Liu, C., Xu, M. M., Zhao, S., Shen, B., Gao, Y., Han, D., He, C. **N6‑methyladenosine of chromosome‑associated regulatory RNA regulates chromatin state and transcription.** *Science* (2020). https://pubmed.ncbi.nlm.nih.gov/31949099/

FYI: I ran it on MacBook Pro 4, there might be differences between macOS, Windows, and Linux.

## Data Download

In Kang et al. (2024), the raw sequencing data were deposited in GEO [GSE230880](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230880) and can be downloaded using SRA accession numbers SRR24339407-SRR24339410.

Sample: (human chromatin-associated A375 cells): 
- METTL3 knockdown (shMETTL3_Input [MeRIP])
    - [SRR24339408](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRX20140110&o=acc_s%3Aa)  shMETTL3_Input_rep1
    - [SRR24339407](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRX20140112&o=acc_s%3Aa)  shMETTL3_Input_rep2
- Control (shControl_Input [MeRIP])
    - [SRR24339410](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRX20140106&o=acc_s%3Aa) shControl_Input_rep1
    - [SRR24339409](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRX20140108&o=acc_s%3Aa)  shControl_Input_rep2

To download the SRA files, install the SRA Toolkit from this github page https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit. For convenience (and to show you where the binaries are) append the path to the binaries to your PATH environment variable: 
```bash
# Example for bash/zsh (adjust path to your installation)
export PATH="$HOME/sratoolkit.3.0.0-mac64/bin:$PATH"
```
Once installed, you can use `prefetch` to download the SRA files and `fasterq-dump` to convert them into FASTQ format.

```bash
#set file path
GEOROOT = "$HOME/Desktop/Li_Lab/GEO_data"
mkdir -p "$GEOROOT/sra_cache" "$GEOROOT/fastq"

# List of accessions
ACCESSIONS=("SRR24339407" "SRR24339408" "SRR24339409" "SRR24339410")

# Loop through each accession
for ACCESSION in "${ACCESSIONS[@]}"; do
    echo "Processing $ACCESSION ..."

    # Download raw SRA file to cache
    prefetch $ACCESSION -O "$GEOROOT/sra_cache"

    # Convert to FASTQ (paired-end split, 8 threads, pipeline mode)
    fasterq-dump $ACCESSION --split-files -e 8 -p -O "$GEOROOT/fastq"

    # Compress FASTQ files
    gzip "$GEOROOT/fastq/${ACCESSION}"_*.fastq

    # Check downloaded files
    ls -lh "$GEOROOT/fastq" | grep $ACCESSION
done
```
## Pipeline
### 1. FastQC check before trimming
Run FastQC on the raw paired FASTQ files to assess base quality, adapter content, overrepresented sequences, etc. This gives you a baseline to compare against post-trim QC.

- Install FastQC [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). It requires a Java environment.

```bash
# Paths (adjust once)
GEOROOT="$HOME/Desktop/Li_Lab/GEO_data/fastq"
FASTQC="$HOME/Desktop/FastQC/fastqc"   # FastQC executable

# Choose a sample
ACCESSION=SRR24339407

# Make output folder for reports
mkdir -p "$GEOROOT/${ACCESSION}/fastqc_pre"

# Run FastQC on raw FASTQ (paired)
"$FASTQC" -t 8 \
  "$GEOROOT/${ACCESSION}/${ACCESSION}_1.fastq.gz" \
  "$GEOROOT/${ACCESSION}/${ACCESSION}_2.fastq.gz" \
  -o "$GEOROOT/${ACCESSION}/fastqc_pre"
```
 Open the generated *_fastqc.html files to review per-base quality and adapter warnings before trimming.

### 2. Trimming data using Trimmomatic
Trim adapters and low-quality bases using Trimmomatic (paired-end mode). We use TruSeq3 adapter set and modest quality filters commonly used for RNA-seq. Use **paired-end** sequence for further analysis.
- Install Trimmomatic [here](https://github.com/usadellab/Trimmomatic)
```bash
# Paths
GEOROOT="$HOME/Desktop/Li_Lab/GEO_data/fastq"
TRIMJAR="$HOME/Desktop/Trimmomatic-0.40/trimmomatic-0.40.jar"
ADAPTERS="$HOME/Desktop/Trimmomatic-0.40/adapters"

# Sample
ACCESSION=SRR24339407
cd "$GEOROOT/${ACCESSION}"

# Run Trimmomatic (paired-end)
# - ILLUMINACLIP removes adapters (TruSeq3-PE.fa)
# - LEADING/TRAILING remove low-Q ends
# - SLIDINGWINDOW trims when avg Q < 15 in 4-bp window
# - MINLEN drops reads <36bp
java -jar "$TRIMJAR" PE -threads 8 \
  "${ACCESSION}_1.fastq.gz" "${ACCESSION}_2.fastq.gz" \
  "${ACCESSION}_1_paired.fq.gz"   "${ACCESSION}_1_unpaired.fq.gz" \
  "${ACCESSION}_2_paired.fq.gz"   "${ACCESSION}_2_unpaired.fq.gz" \
  ILLUMINACLIP:"$ADAPTERS/TruSeq3-PE.fa":2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
  -summary "${ACCESSION}_trim_summary.txt"

# Organize trimmed outputs
mkdir -p "${ACCESSION}_trim_data"
mv "${ACCESSION}_1_paired.fq.gz"  "${ACCESSION}_2_paired.fq.gz" \
   "${ACCESSION}_1_unpaired.fq.gz" "${ACCESSION}_2_unpaired.fq.gz" \
   "${ACCESSION}_trim_summary.txt" "${ACCESSION}_trim_data/"

# Quick check
ls -lh "${ACCESSION}_trim_data"
```

### 3. FastQC check after trimming
Re-run FastQC on the **paired trimmed reads** to confirm that adapter warnings are gone and base-quality profiles improved.

```bash
GEOROOT="$HOME/Desktop/Li_Lab/GEO_data/fastq"
FASTQC="$HOME/Desktop/FastQC/fastqc"

ACCESSION=SRR24339407

# Output folder for post-trim QC
mkdir -p "$GEOROOT/${ACCESSION}/fastqc_post"

# Run FastQC on trimmed paired files only
"$FASTQC" -t 8 \
  "$GEOROOT/${ACCESSION}/${ACCESSION}_trim_data/${ACCESSION}_1_paired.fq.gz" \
  "$GEOROOT/${ACCESSION}/${ACCESSION}_trim_data/${ACCESSION}_2_paired.fq.gz" \
  -o "$GEOROOT/${ACCESSION}/fastqc_post"
```

### 4. Alignment with HISAT2
 Align trimmed paired reads to the hg38 genome using HISAT2 (spliced aligner). Libraries are reverse-stranded, so we use --rna-strandness RF. Output is streamed to samtools sort to produce a coordinate-sorted BAM, then we index and record alignment stats.





