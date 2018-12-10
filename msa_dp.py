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
    python msa_dp.py -f sequences.fasta -s score_matrix.json -d 100
'''

import argparse
import json
import numpy as np
import itertools

'''
Returns tuple representing difference of t1 and t2.
Assumes they are the same length.
'''
def tuple_diff(t1,t2):
    ret = []
    for i in range(len(t1)):
        ret.append(t1[i]-t2[i])
    return tuple(ret)

'''Computes the actual string alignments given the traceback matrix.
Arguments:
    sequences: list of sequences to align
    t: the traceback matrix
Returns:
    aligned: list of aligned sequences
'''
def traceback(sequences, t):
    curr = []
    for s in sequences:
        curr.append(len(s))
    curr = tuple(curr)
    acc = [""] * len(sequences)
    zero = tuple(np.zeros(len(sequences)))
    while curr != zero:
        transition = t[curr]
        next_point = []
        for i in range(len(transition)):
            if transition[i] == 1:
                next_point.append(curr[i]-1)
                acc[i] += sequences[i][curr[i]-1]
            else:
                next_point.append(curr[i])
                acc[i] += '-'
        curr = tuple(next_point)
    ret = []
    for s in acc:
        ret.append(''.join(reversed(s)))
    return ret
    
'''
Generates list of valid transitions in the matrix that could lead to p.
Filters out impossible transitions - out of bounds i.e.
'''
def generate_valid_befores(p,n):
    ret = []
    lst = list(itertools.product([0, 1], repeat=n))
    for t in lst:
        flag = True
        for i in range(len(t)):
            if p[i]-t[i] < 0:
                flag = False
        if flag:
            ret.append(t)
    ret.remove(tuple(np.zeros(n)))
    return ret

'''
Gets the bases for the corresponding transition and point tuples.
Returns a list of bases, and '-' if gap
'''
def get_bases(sequences, transition, point):
    ret = []
    for i in range(len(transition)):
        if transition[i] == 1:
            ret.append(sequences[i][point[i]-1])
        else:
            ret.append('-')
    return ret

'''
Calculates score of column in multiple alignment.
'''
def get_score(col, s, d):
    score = 0
    combinations = itertools.combinations(col,2)
    for c in combinations:
        if c[0] == '-' and c[1] == '-':
            score += 0
        elif c[0] == '-' or c[1] == '-':
            score -= d
        else:
            score += s[c[0]][c[1]]
    return score 

'''
Computes the score and alignment of sequences.
Arguments:
    sequences: list of sequences to align
    s: the score matrix
    d: the gap opening/extension penalty
Returns:
    score: the score of the optimal sequence alignment
    aligned: list of aligned sequences
'''
def sequence_alignment(sequences, scores, d):
    # initialize traceback and dp matrices
    n = len(sequences)
    normal_lengths = []
    lengths = []
    for s in sequences:
        normal_lengths.append(len(s))
        lengths.append(len(s)+1)
    t = np.empty(lengths, tuple)
    dp = np.empty(lengths, np.int32)

    # create sorted list of points to iterate through
    ranges = []
    for l in lengths:
        ranges.append(list(range(l)))
    pointIterator = itertools.product(*ranges)

    # recursive step
    for point in pointIterator:
        # base case
        if point == tuple(np.zeros(n)):
            dp[point] = 0
            continue
        # generate list of transitions that lead to this point
        befores = generate_valid_befores(point, n)
        # find max
        max_transition = befores[0]
        before_point = tuple_diff(point, max_transition)
        current_bases = get_bases(sequences, max_transition, point)
        max_score = dp[before_point] + get_score(current_bases,scores,d)
        for transition in befores:
            before_point = tuple_diff(point, transition)
            current_bases = get_bases(sequences, transition, point)
            score = dp[before_point] + get_score(current_bases,scores,d)
            if score > max_score:
                max_transition = transition
                max_score = score 
        # update dp and traceback
        dp[point] = max_score 
        t[point] = max_transition
        
    aligned = traceback(sequences, t)
    #return max_score, and list of alignments produced by algorithm
    return dp[tuple(normal_lengths)], aligned
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
   
    score, aligned = sequence_alignment(sequences, s, d)
    print "Alignment:"
    print_alignment(aligned)
    print "Score: " + str(score)

if __name__ == "__main__":
    main()
