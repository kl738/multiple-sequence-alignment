'''
Provides helper functions is using heuristics.
'''

import numpy as np

# sequence weighting

# initial gap penalties
def correctGOP(d, s, a1, a2):
    # calculate avg non-diagonal values of score matrix
    res_mismatch = 0
    count = 0
    for i in s:
        for j in s[i]:
            if i != j:
                res_mismatch += s[i][j]
                count += 1
    res_mismatch = res_mismatch/float(count)
    # calculate lengths
    n = len(a1[0])
    m = len(a2[0])
    min_length = min(n,m)
    # calculate percent identity
    total = 0
    count = 0
    for i in range(min_length):
        hasGap = False
        allEqual = True
        compare = a1[0][i]
        for j in range(1,len(a1)):
            if a1[j][i] == '-':
                hasGap = True 
            if a1[j][i] != compare:
                allEqual = False
        for j in range(len(a2)):
            if a2[j][i] == '-':
                hasGap = True 
            if a2[j][i] != compare:
                allEqual = False
        if hasGap:
            continue 
        if allEqual:
            total += 1
            count += 1
        else:
            count += 1
    identity = total/float(count)
    return (d+np.log(min_length)) * res_mismatch * identity

def correctGEP(e, a1, a2):
    log_val = np.log(len(a1[0])/len(a2[0]))
    return (e * (1.0 + abs(log_val)))

# position specific gap penalties
def decreaseOnGap(a, d, e):
    n = len(a[0])
    newD = []
    newE = []
    for i in range(n):
        countGaps = 0
        count = 0
        for s in a:
            if s[i] == '-':
                countGaps += 1
                count += 1
            else:
                count += 1
        if countGaps > 0:
            newD.append(d*0.3*(count-countGaps)/float(count))
            newE.append(e/2.0)
        else:
            newD.append(d)
            newE.append(e)
    return newD, newE
def increaseNearGap(a, oldD):
    n = len(a[0])
    newD = [0] * n
    posWithGap = set()
    for i in range(n):
        for s in a:
            if s[i] == '-':
                posWithGap.add(i)
    posWithoutGap = set(range(n)) - posWithGap
    diffs = [-1,1,-2,2,-3,3,-4,4,-5,5,-6,6,-7,7,-8,8]
    for pos in range(n)
        if pos in posWithoutGap:
            for diff in diffs:
                if pos + diff in posWithGap:
                    newD[pos] = (oldD[pos]*(2+((8-abs(diff))*2)/8))
                    break
        else:
            newD[pos] = oldD[pos]
    return newD
    
# weight matrices

# divergent sequences: Assume done by neighbor-joining already. 

