# fuzzyfastq-rs

## Introduction
The fuzzyfastq is a command-line tool written in Rust to processes FASTQ files (including gzipped FASTQ files) to identify and count reads that match given nucleotide sequences. This tool supports mismatch tolerance, allowing users to specify the percentage of mismatches allowed in the sequence matching process. Results are reported as raw counts and a percentage of total reads. It is intended for determining presence of sequence components and does not take into account multiple sequence matches on a single read.

## Features
Process standard and gzipped FASTQ files.
Match sequences with an allowance for mismatches.
Handle both direct sequence input and sequences provided in a CSV file.
Efficient processing suitable for bioinformatics data analysis.

## Installation

`git clone https://github.com/MLKaufman/fastq-rs`  

`cd fastq-rs`

`cargo install --path .`

## Usage

`fuzzyfastq <mode> <sequence_or_path_to_csv> <fastq_directory> [mismatch_percentage]`

The tool accepts the following command line arguments:

Mode (--seq or --csv): Specifies the mode of operation.
--seq: Directly use a provided sequence.
--csv: Use sequences from a specified CSV file.  

Sequence or Path to CSV File: Depending on the mode, provide either a nucleotide sequence or the path to a CSV file containing sequences.  

FASTQ Directory: Path to the directory containing FASTQ files. 
Input as FASTQ files (.fastq, .fq, .fastq.gz, or .fq.gz formats). 

Mismatch Percentage (optional): The allowable mismatch percentage as a decimal (e.g., 0.1 for 10% mismatches). Defaults to 0 if not provided.


### CSV format
`#Name, Sequence`
`Barcode01,ATGCTACGCTAGCTACGTCAGTCGAT`
`Barcode02,TGCTCGCTAGTCGCATCGATCGATCG`