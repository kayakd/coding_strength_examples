# GEO FASTQ Downloader

This project provides a robust Bash script that downloads FASTQ files for any GEO Series accession (GSE ID) by traversing the full GSE → GSM → SRX → SRR hierarchy via NCBI Entrez Direct tools and the SRA Toolkit.

## Requirements

This script requires:

- **Entrez Direct** (`esearch`, `elink`, `efetch`, `xtract`)  
  → Used to traverse GSE → GSM → SRR via NCBI
- **SRA Toolkit** (`prefetch`, `fasterq-dump`)  
  → Used to download and extract FASTQ files

## Installation Guide

### Recommended: Install via Conda (all platforms)

```bash
conda create -n geo_download_env -c bioconda entrez-direct sra-tools -y
conda activate geo_download_env
```

> This will install all required tools in an isolated environment, with no system-wide impact.

### Linux (Manual Installation)

#### Entrez Direct:
```bash
sh -c "$(wget -q -O - https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
echo 'export PATH="$HOME/edirect:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

#### SRA Toolkit:
```bash
conda install -c bioconda sra-tools
```

### macOS

Same as Linux. Works best with Conda:
```bash
conda install -c bioconda entrez-direct sra-tools
```

> For manual installs: ensure Homebrew is installed and use `brew install sratoolkit` if available.

### Windows

#### Option 1: Use Windows Subsystem for Linux (WSL)

1. [Install WSL](https://learn.microsoft.com/en-us/windows/wsl/install)
2. Open Ubuntu or Debian terminal
3. Then follow **Linux instructions** above using `conda`

#### Option 2: Use Anaconda Prompt (less recommended)

```bash
conda install -c bioconda entrez-direct sra-tools
```

> **Important**: `prefetch` and `fasterq-dump` are command-line tools; run from Git Bash or WSL for full compatibility.

## How to Run the Script

1. Make it executable:
```bash
chmod +x download_gse_fastq.sh
```

2. Run with a GEO Series ID:
```bash
./download_gse_fastq.sh GSE123456
```

## Output

- FASTQ files go to:
  ```
  ~/GEO_DATA/GSE123456/
  ```

- Metadata is saved to:
  ```
  ~/GEO_DATA/GSE123456_RunInfo.csv
  ```
