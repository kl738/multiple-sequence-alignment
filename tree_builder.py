import numpy as np 
'''Computes the actual string alignments given the traceback matrix.
Arguments:
    x: the first string we're aligning
    y: the second string we're aligning
    t: the traceback matrix, which stores values that point to which
       prior matrix was used to reach a given location in each of the
       3 matrices.
    start: value indicating the starting matrix (that had the optimal value)
Returns:
    a_x: the string for the alignment of x's sequence
    a_y: the string for the alignment of y's sequence
'''
def traceback(x, y, t, start):
    a_x = []
    a_y = []
    i = len(x)
    j = len(y)
    while i > 0 or j > 0:
        direction, table = t[i][j][start]
        if direction == 'g':
            a_x.append(x[i-1])
            a_y.append(y[j-1])
            i -= 1
            j -= 1
        elif direction == 'l':
            a_x.append('-')
            a_y.append(y[j-1])
            j -= 1
        elif direction == 'u':
            a_x.append(x[i-1])
            a_y.append('-')
            i -= 1
        start = table
    return (''.join(reversed(a_x)),''.join(reversed(a_y)))

'''
Pairwise sequence alignment from HW2 to determine scores and thus distances.
Used in computing distance matrix.
'''
def compute_distance(x, y, s, d, e):
    m = [[0] * (len(y)+1) for _ in range(len(x)+1)]
    i_x = [[0] * (len(y)+1) for _ in range(len(x)+1)]
    i_y = [[0] * (len(y)+1) for _ in range(len(x)+1)]
    t = []
    for i in range(len(x)+1):
        temp = []
        for j in range(len(y)+1):
            temp.append([['',0],['',0],['',0]])
        t.append(temp)
    for i in range(1,len(x)+1):
        m[i][0] = -1e10
        i_x[i][0] = i_x[i-1][0]-e
        i_y[i][0] = -1e10
        t[i][0][0] = ['u', 1]
        t[i][0][1] = ['u', 1]
        t[i][0][2] = ['u', 1]
    for j in range(1,len(y)+1):
        m[0][j] = -1e10
        i_x[0][j] = -1e10
        i_y[0][j] = i_y[0][j-1]-e
        t[0][j][0] = ['l', 2]
        t[0][j][1] = ['l', 2]
        t[0][j][2] = ['l', 2]
    for i in range(1,len(x)+1):
        for j in range(1,len(y)+1):
            m_m = m[i-1][j-1] + s[x[i-1]][y[j-1]]
            m_ix = i_x[i-1][j-1] + s[x[i-1]][y[j-1]]
            m_iy = i_y[i-1][j-1] + s[x[i-1]][y[j-1]]
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
    a_x, a_y = traceback(x, y, t, start)
    count = 0
    identity = 0
    for i in range(len(a_x)):
        if a_x[i] == '-' or a_y[i] == '-':
            continue
        elif a_x[i] == a_y[i]:
            identity += 1
            count += 1
        else:
            count += 1
    return identity/float(count)


'''
Computes all pairwise distances between sequences using Needlemen-Wunsch.
Returns a distance matrix in the form of dict of dicts.
'''
def compute_distances(sequences,s,d,e):
    n = len(sequences)
    m = np.ndarray((n,n),np.float)
    for i in range(n):
        for j in range(n):
            m[i][j] = compute_distance(sequences[i], sequences[j],s,d,e)
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