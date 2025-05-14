#!/bin/bash

# -----------------------------------------------------------------------------
# Script: download_gse_fastq.sh
# Created by: Koray Dogan Kaya
#
# Description:
# Given a GEO Series accession (e.g., GSE123456), this script:
#   1. Downloads the corresponding SRA RunInfo metadata.
#   2. Extracts GSM sample titles and dynamically fetches SRR run IDs.
#   3. Downloads FASTQ files using prefetch and fasterq-dump (SRA Toolkit).
#   4. Organizes outputs by LIBRARY_STRATEGY (e.g., RNA-Seq).
#   5. Automatically renames files based on GSM title.
#
# Requirements:
#   - Entrez Direct: esearch, efetch, elink, xtract
#   - SRA Toolkit: prefetch, fasterq-dump
#   - Internet access
#   - Core UNIX tools: wget, grep, awk, tr
#
# Usage:
#   ./download_gse_fastq.sh GSE123456
# -----------------------------------------------------------------------------

set -e

GSE_ID=$1
if [ -z "$GSE_ID" ]; then
  echo "Usage: $0 GSE_ID (e.g., GSE157129)"
  exit 1
fi

BASE_DIR=~/GEO_DATA
mkdir -p "${BASE_DIR}/${GSE_ID}"
cd "${BASE_DIR}/${GSE_ID}"

# Fetch GSE docsum XML (contains all GSM sample metadata and their titles/types)
echo "[INFO] Fetching GSE metadata..."
DOCSUM=$(efetch -db gds -id "$GSE_ID" -format docsum | sed -n '/<Samples>/,/<\/Samples>/p')

# Extract all GSM blocks
mapfile -t GSM_BLOCKS < <(echo $DOCSUM | grep -oP '<Sample>.*?</Sample>')

# Fetch SRR mapping per GSM
declare -A GSM_META

for BLOCK in "${GSM_BLOCKS[@]}"; do
  GSM=$(echo "$BLOCK" | grep -oP '(?<=<Accession>)GSM[0-9]+(?=</Accession>)')
  SRR=$(esearch -db sra -query "$GSM" | efetch -format docsum | xtract -pattern DocumentSummary -element Run@acc)
  [[ -z "$SRR" ]] && SRR="NA"

  TITLE=$(echo "$BLOCK" | grep -oP '(?<=<Title>)[^<]+')
  TAG=$(esearch -db sra -query "$GSM" | efetch -format docsum | xtract -pattern DocumentSummary -element LIBRARY_STRATEGY)
  [[ -z "$TAG" ]] && TAG="Unknown"
  TAG=$(echo "$TAG" | tr ' ' '_' | tr -dc '[:alnum:]_-')

  mkdir -p "$TAG"
  cd "$TAG"

  GSM_META[$GSM]="$SRR|$TITLE"
  echo "$GSM|$TITLE|$SRR|$TAG"

  if [[ "$SRR" != "NA" ]]; then
    echo "[INFO] Downloading $SRR for $GSM"
    prefetch --max-size 100G "$SRR"
    fasterq-dump "$SRR" -O . --split-files
    rm -rf "$SRR"

    SAFE_NAME=$(echo "$TITLE" | tr ' ' '_' | tr -dc '[:alnum:]_')
    if [[ -f "${SRR}_1.fastq" && -f "${SRR}_2.fastq" ]]; then
      mv "${SRR}_1.fastq" "${SAFE_NAME}_R1.fastq"
      mv "${SRR}_2.fastq" "${SAFE_NAME}_R2.fastq"
    elif [[ -f "${SRR}.fastq" ]]; then
      mv "${SRR}.fastq" "${SAFE_NAME}.fastq"
    elif [[ -f "${SRR}_1.fastq" ]]; then
      mv "${SRR}_1.fastq" "${SAFE_NAME}.fastq"
    else
      echo "[WARN] No FASTQ output found for $SRR (missing expected FASTQ files)"
      cd ..
    fi
  fi
done
