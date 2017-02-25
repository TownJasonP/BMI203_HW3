from io_sw import *
from smithwaterman import align
import os
import numpy as np

def get_scores(open_pen, cont_pen, matrix):
    '''
    return the positive and negative scores for a given penalty/matrix combo
    '''
    wd = '/Users/student/Desktop/BMIHW3/'
    positive_pairs = read_pairs(wd + 'Pospairs.txt')
    negative_pairs = read_pairs(wd + 'Negpairs.txt')

    sub_code, sub_matrix = matrix

    positive_scores = []
    for i in positive_pairs:
        a = import_fasta(wd + i[0])
        b = import_fasta(wd + i[1])
        score = align(a, b, sub_code, sub_matrix, open_pen, cont_pen)[2]
        positive_scores.append(score)

    negative_scores = []
    for i in negative_pairs:
        a = import_fasta(wd + i[0])
        b = import_fasta(wd + i[1])
        score = align(a, b, sub_code, sub_matrix, open_pen, cont_pen)[2]
        negative_scores.append(score)

    return(positive_scores, negative_scores)

def score_alignments(a, b, matrix):
    num_code, submat = matrix
    assert(len(a) == len(b))
    scores = np.zeros(len(a))
    for i in range(0,len(a)):
        letter_a = num_code[a[i]]
        letter_b = num_code[b[i]]
        scores[i] = submat[letter_a][letter_b] # match score
    return(np.sum(scores))



def ROC(open_pen, cont_pen, matrix, normed):

    wd = os.getcwd() + '/'
    positive_pairs = read_pairs(wd + 'Pospairs.txt')
    negative_pairs = read_pairs(wd + 'Negpairs.txt')

    sub_code, sub_matrix = matrix

    positive_scores = []
    for i in positive_pairs:
        print(i)
        a = import_fasta(wd + i[0])
        b = import_fasta(wd + i[1])
        score = align(a, b, sub_code, sub_matrix, open_pen, cont_pen)[2]
        if normed == True:
            score = score/(min(len(a),len(b)))
        positive_scores.append(score)

    negative_scores = []
    for i in negative_pairs:
        print(i)
        a = import_fasta(wd + i[0])
        b = import_fasta(wd + i[1])
        score = align(a,b, sub_code, sub_matrix, open_pen, cont_pen)[2]
        if normed == True:
            score = score/(min(len(a),len(b)))
        negative_scores.append(score)


    # generate false and true positive pairs from numerous thresholds
    threshholds = np.linspace(0,max(max(positive_scores),max(negative_scores)),10000)
    fp = []
    tp = []

    for t in threshholds:
        true_positives = [i for i in positive_scores if i > t]
        true_positive_rate = len(true_positives)/float(len(positive_scores))
        false_positives = [i for i in negative_scores if i > t]
        false_positive_rate = len(false_positives)/float(len(negative_scores))
        tp.append(true_positive_rate)
        fp.append(false_positive_rate)

    return(fp,tp)

def FPR(open_pen, cont_pen, matrix, TPR):
    '''
    return the FPR for a given gap penalty, matrix combo at a particular TPR
    (TPR between 0 and 1)
    '''
    wd = os.getcwd() + '/'
    positive_pairs = read_pairs(wd + 'Pospairs.txt')
    negative_pairs = read_pairs(wd + 'Negpairs.txt')

    sub_code, sub_matrix = matrix

    positive_scores = []
    for i in positive_pairs:
        a = import_fasta(wd + i[0])
        b = import_fasta(wd + i[1])
        score = align(a, b, sub_code, sub_matrix, open_pen, cont_pen)[2]
        positive_scores.append(score)

    negative_scores = []
    for i in negative_pairs:
        a = import_fasta(wd + i[0])
        b = import_fasta(wd + i[1])
        score = align(a, b, sub_code, sub_matrix, open_pen, cont_pen)[2]
        negative_scores.append(score)

    # calculate false positive rate when true positive rate is 0.7
    # tpr is 0.7 at the 30th percentile of the positive scores dataset
    thresh = np.percentile(positive_scores, 100*(1-TPR))
    false_positives = [i for i in negative_scores if i > thresh]
    fpr = len(false_positives)/float(len(negative_scores))

    return(fpr)

def objective_function(pos, neg):
    '''
    return sum of TP rates for FP rates of 0.0, 0.1, 0.2, and 0.3.
    '''

    fp = []
    tp = []

    for t in [0.0,0.1,0.2,0.3]:
        thresh = np.percentile(neg, 100*(1-t))
        false_positives = [i for i in neg if i > thresh]
        fpr = len(false_positives)/float(len(neg))
        true_positives = [i for i in pos if i > thresh]
        tpr = len(true_positives)/float(len(pos))
        tp.append(tpr)
        fp.append(fpr)

    return(sum(tp))

def ROC_static_align(pos, neg):
        '''
        return TP rates for FP rates within a large range
        '''

        thresholds = np.linspace(0,max(max(pos),max(neg)),10000)
        tp = []
        fp = []
        for thresh in thresholds:
            false_positives = [i for i in neg if i > thresh]
            fpr = len(false_positives)/float(len(neg))
            true_positives = [i for i in pos if i > thresh]
            tpr = len(true_positives)/float(len(pos))
            tp.append(tpr)
            fp.append(fpr)

        return(fp,tp)

def get_alignments(open_pen, cont_pen, matrix):
    '''
    return the positive and negative scores for a given penalty/matrix combo
    '''
    wd = os.getcwd() + '/'
    positive_pairs = read_pairs(wd + 'Pospairs.txt')
    negative_pairs = read_pairs(wd + 'Negpairs.txt')

    sub_code, sub_matrix = matrix

    positive_alignments = []
    for i in positive_pairs:
        a = import_fasta(wd + i[0])
        b = import_fasta(wd + i[1])
        result = align(a, b, sub_code, sub_matrix, open_pen, cont_pen)[0:2]
        positive_alignments.append(result)

    negative_alignments = []
    for i in negative_pairs:
        a = import_fasta(wd + i[0])
        b = import_fasta(wd + i[1])
        result = align(a, b, sub_code, sub_matrix, open_pen, cont_pen)[0:2]
        negative_alignments.append(result)

    return(positive_alignments, negative_alignments)
