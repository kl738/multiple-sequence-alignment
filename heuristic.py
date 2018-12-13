'''
Provides helper functions is using heuristics.
'''

import numpy as np

# sequence weighting
class Node:
    def __init__(self, id = None, left = None, right = None, parent = None, lenToParent = None, num_children = None):
        self.id = id
        self.left = left 
        self.right = right
        self.parent = parent 
        self.lenToParent = lenToParent
        self.num_children = num_children 

def create_node_tree(root, nodes, mapping):
    nodeMap = {}
    def recurse(n):
        if n in mapping:
            return Node(n, num_children=1)
        else:
            lnode,ldist = nodes[n][0]
            rnode,rdist = nodes[n][1]
            left = recurse(lnode)
            right = recurse(rnode)
            node = Node(n, left, right)
            left.parent = node 
            left.lenToParent = ldist
            right.parent = node
            right.lenToParent = rdist
            node.num_children = left.num_children + right.num_children
            nodeMap[n] = node
            return node
    return recurse(root), nodeMap

def get_sequence_weights(mapping,nodeMap):
    ret = []
    for node_id in mapping:
        weight = 0
        node = nodeMap[node_id]
        while node.parent != None:
            weight += node.lenToParent/float(node.num_children)
            node = node.parent
        ret.append(weight)
    return ret
    
# weight matrices
def get_node_distance(nodeMap,n1,n2):
    n1 = nodeMap[n1]
    n2 = nodeMap[n2]
    count1 = 0
    temp1 = n1
    temp2 = n2
    while temp1.parent != None:
        temp1 = temp1.parent
        count1 += 1
    while temp2.parent != None:
        temp2 = temp2.parent
        count2 += 1
    while count1 > count2:
        n1 = n1.parent 
        count1 -= 1 
    while count2 > count1:
        n2 = n2.parent
        count2 -= 1
    dist = 0
    while n1 != n2:
        dist += n1.lenToParent
        dist += n2.lenToParent
        n1 = n1.parent
        n2 = n2.parent 
    return dist
    
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
    for pos in range(n):
        if pos in posWithoutGap:
            for diff in diffs:
                if pos + diff in posWithGap:
                    newD[pos] = (oldD[pos]*(2+((8-abs(diff))*2)/8))
                    break
        else:
            newD[pos] = oldD[pos]
    return newD


# divergent sequences: Assume done by neighbor-joining already. 

