#!/bin/python3

import sys
import re

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
generate_res_dict(headers, sequence)

