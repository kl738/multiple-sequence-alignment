#!/usr/bin/env python

'''
Alignment and score at the bottom of this comment.

Script for computing sequence alignments using multidimensional dynamic
    programming.
Arguments:
    f - FASTA file with sequences in FASTA format.
    s - JSON with the score matrix for alignment.
    d - The gap penalty for the alignment.

Outputs:
    Prints alignment to console.

Example Usage:
    python msa_clustalw.py -f sequences.fasta -s score_matrix.json -d 100
'''

import argparse
import json
import numpy as np
import itertools

'''
Pairwise sequence alignment from HW2 to determine scores and thus distances.
Used in computing distance matrix.
'''
def compute_distance(x, y, s, d):
     # init matrices
    m = [[0] * (len(y)+1) for _ in range(len(x)+1)]
    # base cases
    for i in range(1,len(x)+1):
        m[i][0] = m[i-1][0]-d
    for j in range(1,len(y)+1):
        m[0][j] = m[0][j-1]-d
    # recursive step
    for i in range(1,len(x)+1):
        for j in range(1,len(y)+1):
            vg = m[i-1][j-1] + s[x[i-1]][y[j-1]]
            vu = m[i-1][j] - d
            vl = m[i][j-1] - d
            m[i][j] = max(vg,vu,vl)
    return m[len(x)][len(y)]

'''
Computes all pairwise distances between sequences using Needlemen-Wunsch.
Returns a distance matrix in the form of dict of dicts.
'''
def compute_distances(sequences,s,d):
    n = len(sequences)
    m = np.ndarray((n,n),np.float)
    for i in range(n):
        for j in range(n):
            m[i][j] = compute_distance(sequences[i], sequences[j],s,d)
    d = {}
    for i in range(n):
        temp = {}
        for j in range(len(m[i])):
            temp[j] = m[i][j]
        d[i] = temp
    return d

'''
Builds tree nodes with neighbor-joining, returns nodes and root
'''
def build_tree(D):
    # nodes is the dict of nodes of the tree, z is the id of the new nodes to create >=n
    nodes = {}
    z = len(D)
    def process(z):
        n = len(D)
        # getting row and column sums
        rowsl = {}
        colsl = {}
        for i in D:
            for j in D[i]:
                if rowsl.get(i)==None:
                    rowsl[i] = D[i][j]
                else:
                    rowsl[i]+=D[i][j]
                if colsl.get(j)==None:
                    colsl[j] = D[i][j]
                else:
                    colsl[j]+=D[i][j]
        minv = 1e9
        argmin_i, argmin_j = 0,0
        dij = {i:{} for i in D}
        # creating D*, and fetching min value
        for i in D:
            for j in D[i]:
                if i == j:
                    dij[i][j] = 0
                else:
                    dij[i][j] = (n-2)*D[i][j] - rowsl[i] - colsl[j]
                    if dij[i][j] < minv:
                        minv = dij[i][j]
                        argmin_i = i 
                        argmin_j = j
        # adding edge to nodes
        left_l = .5*(D[argmin_i][argmin_j]+(rowsl[argmin_i]-colsl[argmin_j])/(n-2))
        right_l = .5*(D[argmin_i][argmin_j]+(colsl[argmin_j]-rowsl[argmin_i])/(n-2))
        nodes[z] = []
        nodes[z].append((argmin_i,left_l))
        nodes[z].append((argmin_j,right_l))
        # inserting new row and column
        temp = {}
        for k in D:
            temp[k] = (D[argmin_i][k]+D[argmin_j][k]-D[argmin_i][argmin_j])/2
            D[k][z] = (D[argmin_i][k]+D[argmin_j][k]-D[argmin_i][argmin_j])/2
        temp[z] = 0
        D[z] = temp
        # deleting rows
        D.pop(argmin_i,None)
        D.pop(argmin_j,None)
        # deleting columns
        for i,v in D.items():
            D[i].pop(argmin_i,None)
            D[i].pop(argmin_j,None)   
    while len(D) > 2:
        process(z)
        z += 1
    # add last two
    i,j = D.keys()
    length = D[i][j]
    nodes[z] = [(i,length/2),(j,length/2)]
    return nodes, z

'''
a1 and a2 are lists of sequences and t is the traceback matrix.
Returns the list of joined sequences that are all aligned together.
'''
def traceback(a1, a2, t):
    a1_length = len(a1)
    a2_length = len(a2)
    a1_acc = [[] for _ in range(a1_length)]
    a2_acc = [[] for _ in range(a2_length)]
    i = len(a1[0])
    j = len(a2[0])
    while i > 0 or j > 0:
        # move diagonal
        if t[i][j] == 'g':
            for k in range(a1_length):
                a1_acc[k].append(a1[k][i-1])
            for l in range(a2_length):
                a2_acc[l].append(a2[l][j-1])
            i -= 1
            j -= 1
        # move left
        elif t[i][j] == 'l':
            for k in range(a1_length):
                a1_acc[k].append('-')
            for l in range(a2_length):
                a2_acc[l].append(a2[l][j-1])
            j -= 1
        # move up
        elif t[i][j] == 'u':
            for k in range(a1_length):
                a1_acc[k].append(a1[k][i-1])
            for l in range(a2_length):
                a2_acc[l].append('-')
            i -= 1
    ret = []
    for a in a1_acc:
        ret.append(''.join(reversed(a)))
    for a in a2_acc:
        ret.append(''.join(reversed(a)))
    return ret 

'''
a1 and a2 are lists of sequences that are in their respective alignments
Uses slightly modified pairwise alignment that average scores of a1 and a2
at each recursive step. If gaps are inserted, then they will be inserted
across the entire column within the alignment.
'''
def align(a1,a2,s,d):
    length_a1 = len(a1[0])
    length_a2 = len(a2[0])
    # init matrices
    m = [[0] * (length_a2+1) for _ in range(length_a1+1)]
    t = [[''] * (length_a2+1) for _ in range(length_a1+1)]
    # base cases
    for i in range(1,length_a1+1):
        m[i][0] = m[i-1][0]-d
        t[i][0] = 'u'
    for j in range(1,length_a2+1):
        m[0][j] = m[0][j-1]-d
        t[0][j] = 'l'
    # recursive step
    for i in range(1,length_a1+1):
        for j in range(1,length_a2+1):
            score = 0
            count = 0
            for x in a1:
                for y in a2:
                    if x[i-1] == '-' or y[j-1] == '-':
                        score += 0
                    else:
                        score += s[x[i-1]][y[j-1]]
                    count += 1
            vg = m[i-1][j-1] + score/float(count)
            vu = m[i-1][j] - d
            vl = m[i][j-1] - d
            m[i][j] = max(vg,vu,vl)
            if m[i][j] == vg:
                t[i][j] = 'g'
            elif m[i][j] == vu:
                t[i][j] = 'u'
            elif m[i][j] == vl:
                t[i][j] = 'l'
    aligned = traceback(a1,a2,t)
    return aligned

'''
Uses postorder traversal to align the nodes
'''
def postorder_align(sequences, nodes, root, mapping, s, d):
    def postorder(n):
        lnode,ldist = nodes[n][0]
        rnode,rdist = nodes[n][1]
        if lnode in mapping and rnode in mapping:
            return align([sequences[lnode]],[sequences[rnode]],s,d)
        elif lnode in mapping and rnode not in mapping:
            return align([sequences[lnode]],postorder(rnode),s,d)
        elif lnode not in mapping and rnode in mapping:
            return align(postorder(lnode), [sequences[rnode]],s,d)
        elif lnode not in mapping and rnode not in mapping:
            return align(postorder(lnode),postorder(rnode),s,d)
    aligned = postorder(root)
    return aligned

'''
Computes the score and alignment of sequences.
Arguments:
    sequences: list of sequences to align
    scores: the score matrix
    d: the gap opening/extension penalty
Returns:
    score: the score of the optimal sequence alignment
    aligned: list of aligned sequences
'''
def sequence_alignment(sequences, scores, d):
    #compute distances
    D = compute_distances(sequences, scores, d)
    #compute guide_tree
    nodes, root = build_tree(D)
    # print(nodes,root)
    #postorder traversal
    mapping = dict(enumerate(sequences))
    aligned = postorder_align(sequences, nodes, root, mapping, scores, d)
    return aligned 

'''
Prints aligned sequences in order
'''
def print_alignment(aligned_sequences):
    for s in aligned_sequences:
        print(s)

def main():
    parser = argparse.ArgumentParser(
        description='Calculate sequence alignments for two sequences with a linear gap penalty.')
    parser.add_argument('-f', action="store", dest="f", type=str, required=True)
    parser.add_argument('-s', action="store", dest="s", type=str, required=True)
    parser.add_argument('-d', action="store", dest="d", type=float, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    score_matrix_file = args.s
    d = args.d

    with open(fasta_file) as f:
        sequences = [line.strip() for line in f.readlines()]

    with open(score_matrix_file) as f:
        s = json.loads(f.readlines()[0])
   
    aligned = sequence_alignment(sequences, s, d)
    print "Alignment:"
    print_alignment(aligned)
    # print "Score: " + str(score)

if __name__ == "__main__":
    main()