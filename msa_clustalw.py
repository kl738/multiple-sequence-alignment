#!/usr/bin/env python

'''
Alignment and score at the bottom of this comment.

Script for computing sequence alignments using CLUSTALW
Arguments:
    f - FASTA file with sequences in FASTA format.
    c - configuration file path
        
Outputs:
    Prints alignment to console.

Example Usage:
    python msa_clustalw.py -f sequences.fasta -c config.txt
'''

import argparse
import json
import numpy as np
import itertools
import configuration
import tree_builder
import heuristic

'''
a1 and a2 are lists of sequences and t is the traceback matrix.
Returns the list of joined sequences that are all aligned together.
'''
def traceback(a1, a2, t, start):
    a1_length = len(a1)
    a2_length = len(a2)
    a1_acc = [[] for _ in range(a1_length)]
    a2_acc = [[] for _ in range(a2_length)]
    i = len(a1[0])
    j = len(a2[0])
    while i > 0 or j > 0:
        direction, table = t[i][j][start]
        # move diagonal
        if direction == 'g':
            for k in range(a1_length):
                a1_acc[k].append(a1[k][i-1])
            for l in range(a2_length):
                a2_acc[l].append(a2[l][j-1])
            i -= 1
            j -= 1
        # move left
        elif direction == 'l':
            for k in range(a1_length):
                a1_acc[k].append('-')
            for l in range(a2_length):
                a2_acc[l].append(a2[l][j-1])
            j -= 1
        # move up
        elif direction == 'u':
            for k in range(a1_length):
                a1_acc[k].append(a1[k][i-1])
            for l in range(a2_length):
                a2_acc[l].append('-')
            i -= 1
        start = table
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
def align(a1,a2,config):
    # weight matrix heuristic
    if config.weight_m:
        s = config.scores['BLOSUM62']
    else:
        s = config.scores['BLOSUM62']
    # initial gap penalty heuristic
    if config.init_gap_penalties:
        d = heuristic.correctGOP(config.d, s, a1, a2)
        e = heuristic.correctGEP(config.e, a1, a2)
    else:
        d = config.d
        e = config.e
    # init lengths
    length_a1 = len(a1[0])
    length_a2 = len(a2[0])
    # init matrices
    m = [[0] * (length_a2+1) for _ in range(length_a1+1)]
    i_x = [[0] * (length_a2+1) for _ in range(length_a1+1)]
    i_y = [[0] * (length_a2+1) for _ in range(length_a1+1)]
    t = []
    for i in range(length_a1+1):
        temp = []
        for j in range(length_a2+1):
            temp.append([['',0],['',0],['',0]])
        t.append(temp)    
    # base cases
    for i in range(1,length_a1+1):
        m[i][0] = m[i-1][0]-d
        i_x[i][0] = i_x[i-1][0]-e
        i_y[i][0] = -1e10
        t[i][0][0] = ['u', 1]
        t[i][0][1] = ['u', 1]
        t[i][0][2] = ['u', 1]
    for j in range(1,length_a2+1):
        m[0][j] = m[0][j-1]-d
        i_x[0][j] = -1e10
        i_y[0][j] = i_y[0][j-1]-e
        t[0][j][0] = ['l', 2]
        t[0][j][1] = ['l', 2]
        t[0][j][2] = ['l', 2]
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
            avg_score = score/float(count)
            m_m = m[i-1][j-1] + avg_score
            m_ix = i_x[i-1][j-1] + avg_score
            m_iy = i_y[i-1][j-1] + avg_score
            m[i][j] = max(m_m,m_ix,m_iy)
            if m[i][j] == m_m:
                t[i][j][0] = ['g', 0]
            elif m[i][j] == m_ix:
                t[i][j][0] = ['g', 1]
            elif m[i][j] == m_iy:
                t[i][j][0] = ['g', 2]
            ix_m = m[i-1][j] - d
            ix_ix = i_x[i-1][j] - e
            i_x[i][j] = max(ix_m,ix_ix)
            if i_x[i][j] == ix_m:
                t[i][j][1] = ['u', 0]
            elif i_x[i][j] == ix_ix:
                t[i][j][1] = ['u', 1]
            iy_m = m[i][j-1] - d
            iy_iy = i_y[i][j-1] - e
            i_y[i][j] = max(iy_m,iy_iy)
            if i_y[i][j] == iy_m:
                t[i][j][2] = ['l', 0]
            elif i_y[i][j] == iy_iy:
                t[i][j][2] = ['l', 2]
    max_score = max(m[len(x)][len(y)],i_x[len(x)][len(y)],i_y[len(x)][len(y)])
    if max_score == m[len(x)][len(y)]:
        start = 0
    elif max_score == i_x[len(x)][len(y)]:
        start = 1
    elif max_score == i_y[len(x)][len(y)]:
        start = 2
    aligned = traceback(a1,a2,t,start)
    return aligned

'''
Uses postorder traversal to align the nodes
'''
def postorder_align(sequences, nodes, root, mapping, config):
    def postorder(n):
        lnode,ldist = nodes[n][0]
        rnode,rdist = nodes[n][1]
        if lnode in mapping and rnode in mapping:
            return align([sequences[lnode]],[sequences[rnode]],config)
        elif lnode in mapping and rnode not in mapping:
            return align([sequences[lnode]],postorder(rnode),config)
        elif lnode not in mapping and rnode in mapping:
            return align(postorder(lnode), [sequences[rnode]],config)
        elif lnode not in mapping and rnode not in mapping:
            return align(postorder(lnode),postorder(rnode),config)
    aligned = postorder(root)
    return aligned

'''
Computes the score and alignment of sequences.
Arguments:
    sequences: list of sequences to align
    config: config
Returns:
    score: the score of the optimal sequence alignment
    aligned: list of aligned sequences
'''
def sequence_alignment(sequences, config):
    #compute distances
    D = tree_builder.compute_distances(sequences, config.scores['BLOSUM62'], config.d, config.e)
    #compute guide_tree
    nodes, root = tree_builder.build_tree(D)
    #postorder traversal
    mapping = dict(enumerate(sequences))
    aligned = postorder_align(sequences, nodes, root, mapping, config)
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
    parser.add_argument('-c', action="store", dest="c", type=str, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    config_file = args.c

    with open(fasta_file) as f:
        sequences = [line.strip() for line in f.readlines()]
    
    config = configuration.parseConfigFile(config_file)
    aligned = sequence_alignment(sequences, config)
    print "Alignment:"
    print_alignment(aligned)

if __name__ == "__main__":
    main()