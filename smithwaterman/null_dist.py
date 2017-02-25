'''
Use smithwaterman alignment algorithm on randomly generated strings to look
at null distribution of scores - non-normalized and normalized to the minimum
length of the string
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from io_sw import *
from smithwaterman import align
from utils import *

# import relevant data
wd = os.getcwd() + '/'
positive_pairs = read_pairs(wd + 'Pospairs.txt')
negative_pairs = read_pairs(wd + 'Negpairs.txt')
BLOSUM50 = scoring_matrix(wd + 'BLOSUM50')

# count appearances of each amino acid to get info about frequency
aminos = BLOSUM50[0].keys
dist_p = np.zeros(len(aminos))
dist_n = np.zeros(len(aminos))
# count amino acid appearances in positive pairs
for i in positive_pairs:
    for aa in import_fasta(wd + i[0]):
        for j in enumerate(aminos):
            if aa == j[1]:
                dist_p[j[0]] += 1
    for aa in import_fasta(wd + i[1]):
        for j in enumerate(aminos):
            if aa == j[1]:
                dist_p[j[0]] += 1
# count amino acid appearances in negative pairs
for i in negative_pairs:
    for aa in import_fasta(wd + i[0]):
        for j in enumerate(aminos):
            if aa == j[1]:
                dist_n[j[0]] += 1
    for aa in import_fasta(wd + i[1]):
        for j in enumerate(aminos):
            if aa == j[1]:
                dist_n[j[0]] += 1

dist_t = dist_p + dist_n

# generate a collection of 100 scores for randomly generated strings of size n_a and n_b
# these strings will have realistic amino acid frequencies

scores = []
for iteration in range(100):
    matrix = BLOSUM50
    n_a = 50
    n_b = 50

    seq_a = np.random.choice(aminos, p = dist_t)
    seq_b = np.random.choice(aminos, p = dist_t)

    align_score = align(seq_a, seq_b, BLOSUM50[0], BLOSUM50[1], 6, 4)
    scores.append(align(seq_a, seq_b, BLOSUM50))

expected_score = np.average(scores)
print(expected_score)
