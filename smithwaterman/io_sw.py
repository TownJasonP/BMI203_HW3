import os
import numpy as np

def scoring_matrix(filepath):
    f = open(filepath)
    lines = []
    for line in f:
        if line[0] != '#':
            lines.append(line)

    elements = lines[0].strip('\n').split()
    num_code = dict(zip(elements,range(len(elements))))

    matrix = [i.strip('\n').split() for i in lines[1:]]

    matrix = np.array(matrix, dtype = int)

    return(num_code,matrix)


def import_fasta(filepath):
    f = open(filepath)
    sequence = ''
    for line in f:
        if line[0] != '>':
            sequence += line.strip('\n').upper()
    return(sequence)

def read_pairs(filepath):
    f = open(filepath)
    pairs = []
    for line in f:
        pairs.append(line.strip('\n').split())
    return(pairs)
