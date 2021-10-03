# Introduction to bioinformatics
# Project 1
# PCR Simulation with Python
# Genome E
# 26245-26472
# protein_id="YP_009724392.1"
#"""
#Rob Rouse
#Yuvraj Subedi
#Jaron Smith
#"""

import os
import numpy
import pcr as pcr

print("Reading the fasta file...")
with open('sequence.fasta') as file:
    comments = file.readline()
    SARS_COV2_genome = file.read()

file.close()

# validate the file
SARS_COV2_genome = SARS_COV2_genome.replace('\n', '')

# extract E gene from 26245:26472
E_gene = SARS_COV2_genome[26245:26472]   # rna sequence

# validate gene extracted
E_gene = E_gene.upper()

#get complementary DNA
cDNA_c = pcr.getComp(E_gene)
E_gene = E_gene[::-1]

#double stranded DNA
DNA = (E_gene, cDNA_c)

# get the primers
forwardPrimer = pcr.getPrimers()[0]
reversePrimer = pcr.getPrimers()[1]

# print the primers
print("Forward Primer: " + forwardPrimer[0])
print("Reverse Primer: " + reversePrimer[0])

# Sequences
# print DNA strands to be copied
print(DNA[0])
print(DNA[1])

# PCR cycle begin...
PCR_products = pcr.PCR(DNA, 50, 10)

# Display the statistics of the PCR cycle
pcr.stats(PCR_products)

