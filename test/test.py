#just focused on testing alignment algorithm

from smithwaterman import io
import pytest
import os

filepath = '/Users/student/Desktop/BMIHW3/BLOSUM50'
matrix = scoring_matrix(filepath)

def test_self_alignment():
    # Check that a sequence aligns with itself:
    # first 50 residues of HIF1-alpha
    seq = 'MEGAGGANDKKKISSERRKEKSRDAARSRRSKESEVFYELAHQLPLPHNV'
    align_seq = align(seq,seq,matrix[0],matrix[1],6,4)
    assert align_seq[0] == align_seq[1]
