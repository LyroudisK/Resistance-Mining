#!/bin/python3

import sys
import gzip
import re
#test quim kostas

# Define a function to unzip gzipped files
def g_unzip(gzip_file):
    """
    Gets a gzipped file as input and unzips it.
    Arguments:path to the gzipped file.
    Returns: path to the unzipped file.
    """
    try:
        # Raise an error if provided file is not .gz
        if not gzip_file.endswith('.gz'):
            raise ValueError(f"Error: '{gzip_file}' is not a .gz file.")
        # Unzip file
        unzipped_file = gzip_file.replace('.gz', '')
        with gzip.open(gzip_file, 'rt') as f_in, open(unzipped_file, 'w') as f_out:
            for line in f_in:
                f_out.write(line)
        return unzipped_file
    except FileNotFoundError:
        print(f"Error: File '{gzip_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)


# Define a function to "convert" FASTQ files to FASTA files
def fastq_to_fasta(fastq_file, output_file_name):
    """
    Takes a FASTQ file and stores only the sequence lines to a FASTA file.
    Arguments: FASTQ file (can work also with .txt that has FASTQ content)
               name for the output file
    Returns: a FASTA file containing only sequences but no headers
    """
    output_fasta_file = output_file_name + ".fsa"
    # Open file to read and file to write
    with open(fastq_file, "r") as fr, open(output_fasta_file, "w") as fw:
        dna_flag = False  # Flag to indicate if DNA sequence lines are being read
        for line in fr:
            line = line.strip()
            if line.startswith("@"):  # Start of a new sequence
                dna_flag = True
            elif line.startswith("+"):  # Start of quality scores, end of DNA sequence reading
                dna_flag = False
            elif dna_flag:  # DNA sequence line
                fw.write(line + "\n")
    return output_fasta_file


# Define a function to merge forward and reverse reads
def paired_end_merge(forward_read, reverse_read):
    """
    Takes two FASTA files and merges them.
    Arguments: two FASTA files
    Returns: merged FASTA file
    """
    merged_fasta = "merged_reads.fsa"
    with open(forward_read, "r") as forw, open(reverse_read, "r") as rev, open(merged_fasta, "w") as merged:
        for line in forw:
            merged.write(line)
        for line in rev:
            merged.write(line)
    return merged_fasta


def openfile(filename):
    headers = []
    sequence = []
    with open(filename, "r") as file:
        sequence_line = ""
        sequence_flag = 0
        for line in file:
            line = line.replace("\n", "")
            if line.startswith(">"):
                line = line.replace(">", "")
                line = line.replace(":", "")
                headers.append(line)
                if sequence_flag == 0:
                    pass
                else:
                    sequence.append(sequence_line)
                    sequence_line = ""
                sequence_flag = 1
            else:
                sequence_line += line
    sequence.append(sequence_line)
    return headers, sequence

def generate_kmers(sequence):
    kmers = []
    k = 19
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if len(kmer) == k:  # Check if kmer length is exactly k
            kmers.append(kmer)
    return kmers


def generate_res_dict(headers, sequence):
    res_dict = {}
    # First we calculate the number of sequences that we have to iterate through all.
    for i in range(len(headers)):
        gene = re.split("\s", headers[i])[0]
        resistance = re.split("\s", headers[i])[1]
        kmers = generate_kmers(sequence[i])
        # Check if the resistance key exists, if not, create it
        if resistance not in res_dict:
            res_dict[resistance] = {}
        # Check if the gene key exists under the resistance key, if not, create it
        if gene not in res_dict[resistance]:
            res_dict[resistance][gene] = {}
        # Add kmers to the gene dictionary with a value of 0
        for kmer in kmers:
            res_dict[resistance][gene][kmer] = 0
    count_kmers = 0
    for resistance in res_dict:
        for gene in res_dict[resistance]:
            count_kmers += len(res_dict[resistance][gene].keys())
    return res_dict


headers, sequence = openfile(sys.argv[1])
print(headers)
generate_res_dict(headers, sequence)


if len(sys.argv) == 2:
    forward_read = sys.argv[1]
    f_g_unzip = g_unzip(forward_read)
    f_fasta = fastq_to_fasta(f_g_unzip, 'forward')
elif len(sys.argv) == 3:
    forward_read = sys.argv[1]
    reverse_read = sys.argv[2]
    f_g_unzip = g_unzip(forward_read)
    r_g_unzip = g_unzip(reverse_read)
    f_fasta = fastq_to_fasta(f_g_unzip, 'forward')
    r_fasta = fastq_to_fasta(r_g_unzip, 'reverse')
    merged_fasta = paired_end_merge(f_fasta, r_fasta)
else:
    print("Please provide either one (single-end) or two (paired-end) .gzip FASTQ files.")
    sys.exit(1)
