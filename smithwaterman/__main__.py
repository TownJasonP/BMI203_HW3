import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from io_sw import *
from smithwaterman import align
from utils import *
import time

#import the relevant data
wd = os.getcwd() + '/'
positive_pairs = read_pairs(wd + 'Pospairs.txt')
negative_pairs = read_pairs(wd + 'Negpairs.txt')
BLOSUM50 = scoring_matrix(wd + 'BLOSUM50')

'''
Use smithwaterman alignment algorithm on randomly generated strings to look
at null distribution of scores - non-normalized and normalized to the minimum
length of the string
'''

# count appearances of each amino acid to get info about frequency
# also collect sequence length data

aminos = BLOSUM50[0].keys()
dist = np.zeros(len(aminos))
seq_lengths_pos_1 = []
seq_lengths_pos_2 = []
seq_lengths_neg_1 = []
seq_lengths_neg_2 = []

# count amino acid appearances in positive pairs
for i in positive_pairs:
    seq_lengths_pos_1.append(len(import_fasta(wd + i[0])))
    for aa in import_fasta(wd + i[0]):
        for j in enumerate(aminos):
            if aa == j[1]:
                dist[j[0]] += 1
    seq_lengths_pos_2.append(len(import_fasta(wd + i[1])))
    for aa in import_fasta(wd + i[1]):
        for j in enumerate(aminos):
            if aa == j[1]:
                dist[j[0]] += 1

# count amino acid appearances in negative pairs
for i in negative_pairs:
    seq_lengths_neg_1.append(len(import_fasta(wd + i[0])))
    for aa in import_fasta(wd + i[0]):
        for j in enumerate(aminos):
            if aa == j[1]:
                dist[j[0]] += 1
    seq_lengths_neg_2.append(len(import_fasta(wd + i[1])))
    for aa in import_fasta(wd + i[1]):
        for j in enumerate(aminos):
            if aa == j[1]:
                dist[j[0]] += 1
dist = dist/(np.sum(dist))

# create scatterplot showing sequence length data for
# context of normalization problem
plt.figure(figsize = (6,6), facecolor = 'white')
plt.scatter(seq_lengths_pos_1, seq_lengths_pos_2, color = 'blue', marker = 'o', alpha = 0.7, label = 'positive pairs')
plt.scatter(seq_lengths_neg_1, seq_lengths_neg_2, color = 'red', marker = 'o', alpha = 0.7, label = 'negative pairs')
plt.xlabel('Length Sequence A')
plt.ylabel('Length Sequence B')
plt.axis([0,1000,0,1000])
plt.grid()
plt.legend(loc = 0)
plt.title('Lengths of Paired Sequences')
plt.savefig('lengthdata.png')

# generate a collection of 100 scores for randomly generated strings of size n_a and n_b
# these strings will have realistic amino acid frequencies
ex_score_matrix = np.zeros((8,8))
norm_score_matrix = np.zeros((8,8))
product_norm_score_matrix = np.zeros((8,8))
test_set = range(50,401,50)
for n_a in range(len(test_set)):
    for n_b in range(len(test_set)):
        scores = []
        for iteration in range(25):
            seq_a = np.random.choice(aminos, p = dist, size = test_set[n_a])
            seq_b = np.random.choice(aminos, p = dist, size = test_set[n_b])

            align_score = align(seq_a, seq_b, BLOSUM50[0], BLOSUM50[1], 6, 4)[2]
            scores.append(align_score)
        print(n_a,n_b)
        expected_score = np.average(scores)
        ex_score_matrix[n_a][n_b] = expected_score
        norm_score_matrix[n_a][n_b] = expected_score/min(test_set[n_a],test_set[n_b])
        product_norm_score_matrix[n_a][n_b] = expected_score/((test_set[n_a]**2 + test_set[n_b]**2)**0.5)
plt.figure(figsize = (6,6), facecolor = 'white')
plt.imshow(ex_score_matrix, cmap = 'magma', interpolation = 'nearest')
plt.xticks(range(len(test_set)),test_set, rotation = 'vertical')
plt.yticks(range(len(test_set)),test_set)
plt.xlabel('Length Sequence A')
plt.ylabel('Length Sequence B')
plt.colorbar(shrink = 0.7, label = 'Expected Score')
plt.title('Raw Alignment Scores\nfor Randomly Generated Strings')
plt.savefig('raw_scores_rand.png')

plt.figure(figsize = (6,6), facecolor = 'white')
plt.imshow(norm_score_matrix, cmap = 'magma', interpolation = 'nearest')
plt.xticks(range(len(test_set)),test_set, rotation = 'vertical')
plt.yticks(range(len(test_set)),test_set)
plt.xlabel('Length Sequence A')
plt.ylabel('Length Sequence B')
plt.colorbar(shrink = 0.7, label = 'Normalized Score')
plt.title('Min Normalized Alignment Scores\nfor Randomly Generated Strings')
plt.savefig('norm_scores_rand.png')

plt.figure(figsize = (6,6), facecolor = 'white')
plt.imshow(product_norm_score_matrix, cmap = 'magma', interpolation = 'nearest')
plt.xticks(range(len(test_set)),test_set, rotation = 'vertical')
plt.yticks(range(len(test_set)),test_set)
plt.xlabel('Length Sequence A')
plt.ylabel('Length Sequence B')
plt.colorbar(shrink = 0.7, label = 'Normalized Score')
plt.title('Product Normalized Alignment Scores\nfor Randomly Generated Strings')
plt.savefig('p_norm_scores_rand.png')

'''
Part 1.A
Consider the false positive rate (proportion of negative pairs with scores that exceed a score
threshold) when the true positive rate (proportion of positive pairs with scores above the
threshold) is 0.7. What's the best false positive rate that you can achieve with varying both gap
opening (from 1 to 20) and extension penalties (from 1 to 5) with the BLOSUM50 matrix? What
is the best gap penalty combination?
'''
# generate a histogram of positive and negative scores for a few penatlies to
# get a general idea of the separation
plt.figure(facecolor = 'white', figsize = (8,6), dpi = 300)
count = 1
n = 4
for i in range(1,n):
    for j in range(1,n):
        pos_scores, neg_scores = get_scores(i, j, BLOSUM50)
        plt.subplot(n-1, n-1, count)
        plt.hist(pos_scores, bins = np.linspace(0,100000,26), color = 'blue', alpha = 0.5)
        plt.hist(neg_scores, bins = np.linspace(0,100000,26), color = 'red', alpha = 0.5)
        plt.title('Distribution of Scores\nOpen Penalty: ' + str(i) + ', Extension Penalty: ' + str(j), size = 6)
        plt.xticks([])
        plt.ylim(0,20)
        plt.yticks([])
        count = count + 1
plt.savefig('histograms.png')

if os.path.isfile(wd+'gapdata.npy'):
    print('found gap data')
    fpr_matrix = np.load(wd + 'gapdata.npy')
else:
    print('calculating optimal gap penalties')
    open_penalties = range(1,21)
    ext_penalties = range(1,6)

    fpr_matrix = np.zeros((20,5))
    for o in open_penalties:
        for e in ext_penalties:
            fpr_matrix[o-1][e-1] = FPR(o, e, BLOSUM50, 0.7)
            print(fpr_matrix)
            print('\n')
    #this takes a while, only want to do it once
    np.save('gapdata',fpr_matrix)

plt.figure(facecolor = 'white', figsize = (4,6), dpi = 300)
plt.imshow(fpr_matrix, interpolation = 'nearest', aspect = 1, cmap = 'magma', vmin = 0, vmax = 0.5)
plt.title('Effect of Gap Penalties\non False Positive Rates')
plt.xlabel('continuing gap penalty')
plt.ylabel('opening gap penalty')
plt.xticks(range(5),range(1,6))
plt.yticks(range(20),range(1,21))
plt.colorbar(shrink = 0.8)
plt.savefig('gaps_fpr.png')

'''
Part 1.B
Using the gap penalties you determined from question 1, which of the provided scoring
matrices performs the best, in terms of false positive rate (at a true positive rate of 0.7)? What
are the performance rates of each of the matrices? Create a Receiver Operator Curve (ROC)
graph which shows the fraction of true positives on the Y axis and the fraction of false positives
on the X axis. Include on the graph data for each of the provided matrices. Please take care to
make your ROC graphs square, with both X and Y axes limited to the range [0:1]. Note, you
can download ROC code from here: http://www.jainlab.org/Public/ucsf-roc.zip. It is not
guaranteed to be bug free but it might save you some time.
'''

# get best gap penalty values from fpr_matrix above
(best_open_pen, best_cont_pen) = np.unravel_index(np.argmin(fpr_matrix), dims = fpr_matrix.shape)

# add one to convert to gap penalty from indices
(best_open_pen, best_cont_pen) = (best_open_pen + 1, best_cont_pen + 1)

# import the other substitution matrices
BLOSUM62 = scoring_matrix(wd + 'BLOSUM62')
MATIO = scoring_matrix(wd + 'MATIO')
PAM100 = scoring_matrix(wd + 'PAM100')
PAM250 = scoring_matrix(wd + 'PAM250')

# set up some lists for plotting
sub_schemes = [BLOSUM50,BLOSUM62,MATIO,PAM100,PAM250]
labels = ['BLOSUM50','BLOSUM62','MATIO','PAM100','PAM250']
cmap = cm.get_cmap(name= 'viridis')
colors = [cmap(i) for i in np.linspace(0,1,len(sub_schemes))][::-1]

# align data in positive and negative pair datasets for each substitution scheme
plt.figure(figsize = (6,6), facecolor = 'white', dpi = 300)
for scheme in enumerate(sub_schemes):
    fp,tp = ROC(best_open_pen,best_cont_pen,scheme[1], False)

    plt.plot(fp,tp, lw = 2, color = colors[scheme[0]], alpha = 0.7, label = labels[scheme[0]])
    plt.plot([0,1],[0,1], color = 'k', ls = ':')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc = 0)
    plt.grid()
    plt.axis([-0.05,1.05,-0.05,1.05])
    plt.savefig('ROC.png')


'''
Part 1.C
How does the performance change if you normalize the Smith-Waterman scores by the
length of the shortest sequence in a pair (i.e. divide the raw score by the min length)? Show
the ROC curves for your best matrix and for the same matrix with normalized scores. Are the
false positive rates better or worse? Why do you think this is so?
'''


#optimal conditions as defined by previous code:
open_pen = 6
cont_pen = 4
matrix = BLOSUM50

#generate ROC data for non-normalized and normalized data
fp,tp = ROC(open_pen, cont_pen, matrix, False)
fp_norm,tp_norm = ROC(open_pen, cont_pen, matrix, True)

plt.figure(figsize = (6,6), facecolor = 'white', dpi = 300)
plt.plot(fp,tp, lw = 2, color = 'blue', alpha = 0.7, label = 'non-normalized')
plt.plot(fp_norm,tp_norm, lw = 2, color = 'red', alpha = 0.7, label = 'normalized')
plt.plot([0,1],[0,1], color = 'k', ls = ':')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc = 0)
plt.axis([-0.05,1.05,-0.05,1.05])
plt.grid()
plt.savefig('NormedScoresROC.png')



'''
Using the best gap penalties and matrix from part 1, create an alignment for each positive pair
of sequences and each negative pair of sequences. You will use these static alignments as a
starting point from which to optimize a scoring matrix to maximize separation of scores of the
positive and negative pairs.
'''

positive_alignments, negative_alignments = get_alignments(6, 4, BLOSUM50)


'''
1. Devise an optimization algorithm to modify the values in a starting score matrix such as to
maximize the following objective function: sum of TP rates for FP rates of 0.0, 0.1, 0.2, and
0.3. The maximum value for the objective function is 4.0 (where you are getting perfect
separation of positive and negative pairs even at the lowest false positive rate). You should
use the gap and extension penalties derived from Part 1. Remember, you must maintain
symmetry in your matrix. You can make use of real-valued scores in the matrices if desired
(this is probably a good idea).
'''

#See utils.py for objective function

'''
2. Beginning from the best matrix from above (that which produced the alignments), run your
optimization algorithm to maximize the fitness of the new matrix. How much improvement do
you see in the fitness? Show the full ROC curves for the original matrix and the optimized
matrix. What happens when you now realign the sequences using the new matrix and
rescore? Show the new ROC curve following realignment on the same graph as above.
Qualitatively discuss how your matrix and your alignments change following optimization.
'''

open_pen = 6
cont_pen = 4
best_mutant = (BLOSUM50[0],BLOSUM50[1].copy())
positive = [score_alignments(i[0], i[1], best_mutant) for i in positive_alignments]
negative = [score_alignments(i[0], i[1], best_mutant) for i in negative_alignments]
obj = objective_function(positive, negative)

#keep track of any improvement in the objective function over each iteration
epoch = [0]
obj_score = [obj]
n = len(best_mutant[1])

generations = 100

for iteration in range(generations):

    start_time = time.time()

    max_fitness = 0

    #generate a population of 'mutated' matrices
    for mutant in range(1000):
        #generate a 'mutated' matrix with mutation_rate differences
        mutation_rate = 3
        pos = np.random.randint(0, n, size = (mutation_rate,2))
        k = np.random.choice([-1,1], size = mutation_rate)
        mutated_matrix = (best_mutant[0], best_mutant[1].copy())


        for i in range(mutation_rate):
            #keep matrix symmetric if values are off diagonal
            if pos[i][0] != pos[i][1]:
                mutated_matrix[1][pos[i][0],pos[i][1]] = mutated_matrix[1][pos[i][0]][pos[i][1]] + k[i]
                mutated_matrix[1][pos[i][1],pos[i][0]] = mutated_matrix[1][pos[i][1]][pos[i][0]] + k[i]
            else:
                mutated_matrix[1][pos[i][0],pos[i][1]] = mutated_matrix[1][pos[i][0]][pos[i][1]] + k[i]

        # calculate objective function for mutated matrix
        positive = [score_alignments(a[0], a[1], mutated_matrix) for a in positive_alignments]
        negative = [score_alignments(a[0], a[1], mutated_matrix) for a in negative_alignments]
        obj = objective_function(positive, negative)
        if obj > max_fitness:
            max_fitness = obj
            best_mutant = (mutated_matrix[0],mutated_matrix[1].copy())

    epoch.append(iteration + 1)
    obj_score.append(max_fitness)
    print('epoch ' + str(iteration+1) + ' runtime: ' + str(time.time()-start_time))
    print('max fitness this round: ' + str(max_fitness))
    print(BLOSUM50[1] - best_mutant[1])
    #select most fit 'mutant' to become new starting matrix


plt.figure(figsize = (6,6), facecolor = 'white', dpi = 300)
plt.plot(epoch,obj_score, color = 'blue', alpha = 0.7)
plt.axhline(4, color = 'black', ls = ':')
plt.xlabel('Epoch')
plt.ylabel('Objective Function Score')
plt.title('Genetic Optimization of the BLOSUM50 Matrix')
plt.axis([0,len(epoch),0,5])
plt.savefig('BLOSUM50_Optimization.png')

#for static alignments
positive = [score_alignments(a[0], a[1], BLOSUM50) for a in positive_alignments]
negative = [score_alignments(a[0], a[1], BLOSUM50) for a in negative_alignments]
x1,y1 = ROC_static_align(positive,negative)

positive = [score_alignments(a[0], a[1], best_mutant) for a in positive_alignments]
negative = [score_alignments(a[0], a[1], best_mutant) for a in negative_alignments]
x2,y2 = ROC_static_align(positive,negative)


plt.figure(figsize = (6,6), facecolor = 'white', dpi = 300)
plt.plot(x1,y1, color = 'blue', alpha = 0.7, lw = 2, label = 'BLOSUM50')
plt.plot(x2,y2, color = 'red', alpha = 0.7, lw = 2, label = 'Optimized BLOSUM50')
plt.plot([-0,1],[-0,1], color = 'k', ls = ':')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Static Alignment with Optimized BLOSUM50 Matrix')
plt.axis([-0.05, 1.05, -0.05, 1.05])
plt.legend(loc = 0)
plt.grid()
plt.savefig('Optimized_BLOSUM50_ROC_Static.png')


# realign using new matrix and evaluate ROC curve
x1,y1 = ROC(open_pen, cont_pen, BLOSUM50, False)
x2,y2 = ROC(open_pen, cont_pen, best_mutant, False)

plt.figure(figsize = (6,6), facecolor = 'white', dpi = 300)
plt.plot(x1,y1, color = 'blue', alpha = 0.7, lw = 2, label = 'BLOSUM50')
plt.plot(x2,y2, color = 'red', alpha = 0.7, lw = 2, label = 'Optimized BLOSUM50')
plt.plot([0,1],[0,1], color = 'k', ls = ':')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc = 0)
plt.axis([-0.05, 1.05, -0.05, 1.05])
plt.title('Nonstatic Alignment with Optimized BLOSUM50 Matrix')
plt.grid()
plt.savefig('Optimized_BLOSUM50_ROC_alignment.png')

'''
3. Beginning from the MATIO matrix, but using the same initial sequence alignments, re-run
the optimization. Show the same ROC plots as for (2). Discuss the relationship between the
results you see here and the results you saw for (2).
'''

open_pen = 6
cont_pen = 4
best_mutant = (MATIO[0],MATIO[1].copy())
positive = [score_alignments(i[0], i[1], best_mutant) for i in positive_alignments]
negative = [score_alignments(i[0], i[1], best_mutant) for i in negative_alignments]
obj = objective_function(positive, negative)

#keep track of any improvement in the objective function over each iteration
epoch = [0]
obj_score = [obj]
n = len(best_mutant[1])

generations = 100

for iteration in range(generations):

    start_time = time.time()

    max_fitness = 0

    #generate a population of 'mutated' matrices
    for mutant in range(1000):
        #generate a 'mutated' matrix with mutation_rate differences
        mutation_rate = 3
        pos = np.random.randint(0, n, size = (mutation_rate,2))
        k = np.random.choice([-1,1], size = mutation_rate)
        mutated_matrix = (best_mutant[0], best_mutant[1].copy())


        for i in range(mutation_rate):
            #keep matrix symmetric if values are off diagonal
            if pos[i][0] != pos[i][1]:
                mutated_matrix[1][pos[i][0],pos[i][1]] = mutated_matrix[1][pos[i][0]][pos[i][1]] + k[i]
                mutated_matrix[1][pos[i][1],pos[i][0]] = mutated_matrix[1][pos[i][1]][pos[i][0]] + k[i]
            else:
                mutated_matrix[1][pos[i][0],pos[i][1]] = mutated_matrix[1][pos[i][0]][pos[i][1]] + k[i]

        # calculate objective function for mutated matrix
        positive = [score_alignments(a[0], a[1], mutated_matrix) for a in positive_alignments]
        negative = [score_alignments(a[0], a[1], mutated_matrix) for a in negative_alignments]
        obj = objective_function(positive, negative)
        if obj > max_fitness:
            max_fitness = obj
            best_mutant = (mutated_matrix[0],mutated_matrix[1].copy())

    epoch.append(iteration + 1)
    obj_score.append(max_fitness)
    print('epoch ' + str(iteration+1) + ' runtime: ' + str(time.time()-start_time))
    print('max fitness this round: ' + str(max_fitness))
    print(MATIO[1] - best_mutant[1])
    #select most fit 'mutant' to become new starting matrix


plt.figure(figsize = (6,6), facecolor = 'white', dpi = 300)
plt.plot(epoch,obj_score, color = 'blue', alpha = 0.7)
plt.axhline(4, color = 'black', ls = ':')
plt.xlabel('Epoch')
plt.ylabel('Objective Function Score')
plt.title('Genetic Optimization of the MATIO Matrix')
plt.axis([0,len(epoch),0,5])
plt.savefig('MATIO_Optimization.png')

#for static alignments
positive = [score_alignments(a[0], a[1], MATIO) for a in positive_alignments]
negative = [score_alignments(a[0], a[1], MATIO) for a in negative_alignments]
x1,y1 = ROC_static_align(positive,negative)

positive = [score_alignments(a[0], a[1], best_mutant) for a in positive_alignments]
negative = [score_alignments(a[0], a[1], best_mutant) for a in negative_alignments]
x2,y2 = ROC_static_align(positive,negative)


plt.figure(figsize = (6,6), facecolor = 'white', dpi = 300)
plt.plot(x1,y1, color = 'blue', alpha = 0.7, lw = 2, label = 'MATIO')
plt.plot(x2,y2, color = 'red', alpha = 0.7, lw = 2, label = 'Optimized MATIO')
plt.plot([-0,1],[-0,1], color = 'k', ls = ':')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Static Alignment with Optimized MATIO Matrix')
plt.axis([-0.05, 1.05, -0.05, 1.05])
plt.legend(loc = 0)
plt.grid()
plt.savefig('Optimized_MATIO_ROC_Static.png')


# realign using new matrix and evaluate ROC curve
x1,y1 = ROC(open_pen, cont_pen, MATIO, False)
x2,y2 = ROC(open_pen, cont_pen, best_mutant, False)

plt.figure(figsize = (6,6), facecolor = 'white', dpi = 300)
plt.plot(x1,y1, color = 'blue', alpha = 0.7, lw = 2, label = 'MATIO')
plt.plot(x2,y2, color = 'red', alpha = 0.7, lw = 2, label = 'Optimized MATIO')
plt.plot([0,1],[0,1], color = 'k', ls = ':')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc = 0)
plt.axis([-0.05, 1.05, -0.05, 1.05])
plt.title('Nonstatic Alignment with Optimized MATIO Matrix')
plt.grid()
plt.savefig('Optimized_MATIO_ROC_alignment.png')
