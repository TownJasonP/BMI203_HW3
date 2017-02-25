import numpy as np


def align(sequence_a, sequence_b, num_code, submat, open_pen, cont_pen):

    n_a = len(sequence_a)
    n_b = len(sequence_b)

    # Initialize the score matrix
    matrix_d = np.zeros((n_a+1,n_b+1)) #matrix for diagonal scores
    matrix_i = np.zeros((n_a+1,n_b+1)) #matrix for diagonal scores
    matrix_j = np.zeros((n_a+1,n_b+1)) #matrix for diagonal scores
    direction = np.zeros((n_a+1,n_b+1)) #matrix containing directional info

    max_val = 0
    max_pos = None

    for i in range(1,n_a+1):
        for j in range(1,n_b+1):
            letter_a = num_code[sequence_a[i-1]]
            letter_b = num_code[sequence_b[j-1]]

            #recursion relations for matrix_d:
            s = submat[letter_a][letter_b] # match score
            a = matrix_d[i-1][j-1] + s
            b = matrix_i[i-1][j-1] + s
            c = matrix_j[i-1][j-1] + s
            matrix_d[i][j] = max(a,b,c,0)

            direction[i][j] = [0,a,b,c].index(max(a,b,c,0))

            #recursion relations for matrix_i:
            a = matrix_d[i-1][j] - open_pen
            b = matrix_i[i-1][j] - cont_pen
            matrix_i[i][j] = max(a,b)

            #recursion relations for matrix_j:
            a = matrix_d[i][j-1] - open_pen
            b = matrix_i[i][j-1] - cont_pen
            matrix_j[i][j] = max(a,b)

            #record maximum value for traceback starting location
            if matrix_d[i][j] > max_val:
                max_val = matrix_d[i][j]
                max_pos = [i,j]

    # Traceback
    location = max_pos
    alignment_a = sequence_a[location[0]-1] # first letter in alignment
    alignment_b = sequence_b[location[1]-1] # first letter in alignment
    where_to = direction[location[0]][location[1]] # first to direction to follow

    score = 0 # initialize score to zero

    while where_to != 0: # end when we get to a stopping point

        score += matrix_d[location[0]][location[1]]

        if where_to == 1:
            # move diagonally
            location[0] = location[0] - 1
            location[1] = location[1] - 1

            alignment_a += sequence_a[location[0]-1]
            alignment_b += sequence_b[location[1]-1]

            where_to = direction[location[0]][location[1]]

        elif where_to == 3:
            # move left, meaning a gap in protein a
            location[0] = location[0]
            location[1] = location[1] - 1

            alignment_a += '*'
            alignment_b += sequence_b[location[1]-1]

            where_to = direction[location[0]][location[1]]

        elif where_to == 2:
            # move up, meaning a gap in protein b
            location[0] = location[0] - 1
            location[1] = location[1]

            alignment_a += sequence_a[location[0]-1]
            alignment_b += '*'

            where_to = direction[location[0]][location[1]]

    return(alignment_a[::-1], alignment_b[::-1], score)
