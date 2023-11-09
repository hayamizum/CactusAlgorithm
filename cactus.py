#Constructing Cactus
#
#written by Keita Watanabe, Waseda university
#2023.11.9
#
#
import sys
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import itertools
import copy
import csv

#Function readcsv
#Input data from csv file
def readcsv(filename):
    distance_matrix = []
    with open(filename, 'r', newline='', encoding='utf-8-sig') as file:
        reader = csv.reader(file)
        for i, row in enumerate(reader):
            if i == 0:
                #First row is size of matrix
                n = int(row[0])
                continue
            #Add the data from the second line to a distance matrix
            distance_matrix.append([float(value) for value in row])
    return n, np.array(distance_matrix)

#Function comp
#Calculate compactification index in V
def comp(V, x):

    global D
    
    #Calculate minimum value
    if x==list(V)[0]:
        y,z = list(V)[1], list(V)[2]
        min_val = (D[x,y]+D[x,z]-D[y,z])/2
        min_comb = (y,z)
    elif x==list(V)[1]:
        y,z = list(V)[0], list(V)[2]
        min_val = (D[x,y]+D[x,z]-D[y,z])/2
        min_comb = (y,z)
    else:
        y,z = list(V)[0], list(V)[1]
        min_val = (D[x,y]+D[x,z]-D[y,z])/2
        min_comb = (y,z)
    
    for y, z in itertools.combinations(V-{x},2):
        val = (D[x,y]+D[x,z]-D[y,z])/2
        if val < min_val:
            min_val = val
            min_comb = (y, z)
    
    return min_val, min_comb    

##Function Compactification (Tightening)
#Determine how much each vertex "extends" from the set of other vertices using tightening (Compactification).
def Compactification(V):
    global D
    global VertexNum
    global V_del
    global AllVertices
    global Adj
    global Adj_opt
    V_aux = set()

    print(V)
    if len(V)<=2:
        mkAdj(V)
        return
    
    VertexNum_tmp = VertexNum
    doEnd = 0

    for x in V:
        index, (y, z) = comp(V,x)
        if index == 0:
            continue
        else:
            #Auxiliary vertex
            VertexNum += 1
            aux = VertexNum
            #Expand distance matrix
            D = np.append(D, np.zeros((1, D.shape[0])), axis=0)
            D = np.append(D, np.zeros((D.shape[0], 1)), axis=1)
            #Write x, y, z into the distance matrix
            D[x,aux] = D[aux,x] = index
            D[y,aux] = D[aux,y] = D[y,x] - D[x,aux]
            D[z,aux] = D[aux,z] = D[z,x] - D[x,aux]
            #Identification
            isEqualed = 0
            for a in {x,y,z}:
                if D[a, aux] == 0:
                    if a in V_del:
                        doEnd = 1
                    VertexNum -= 1
                    D = np.delete(D, aux, axis=0)
                    D = np.delete(D, aux, axis=1)
                    isEqualed = 1
                    break
            if isEqualed == 1:
                continue
            #Write other verticies into the distance matrix
            for a in AllVertices.union(V_aux)-{x,y,z,aux}:
                D[a,aux] = D[aux,a] = max(D[a,t] - D[t,aux] for t in {x,y,z})
                #Identification
                if D[a, aux] == 0:
                    if a in V_del:
                        doEnd = 1
                    VertexNum -= 1
                    D = np.delete(D, aux, axis=0)
                    D = np.delete(D, aux, axis=1)
                    isEqualed = 1
                    break
            if isEqualed == 0:
                V_aux.add(aux)
                AllVertices.add(aux)

    if VertexNum_tmp==VertexNum and doEnd==1:
        for x,y in itertools.permutations(V,2):
            if Adj[x,y]==1:
                Adj[x,y] = 0
                Adj_opt[x,y] = 1
        V_del = V_del.union(V)
        return
    
    mkAdj(V.union(V_aux))



#Function CountNeighbor
#Count neighbor of x
def CountNeighbor(x):
    global Adj
    NeighborNum = 0
    for i in range(len(Adj)):
        if Adj[x,i]==1:
            NeighborNum += 1
    return NeighborNum


#Function CountOptNeighbor
#Count neighbor of x in optimal realization
def CountOptNeighbor(x):
    global Adj_opt
    NeighborNum = 0
    for i in range(len(Adj_opt)):
        if Adj_opt[x,i]==1:
            NeighborNum += 1
    return NeighborNum


#Function makeAdjacent
def mkAdj(V):
    global D
    global V_del
    global AllVertices
    global Adj
    global Adj_opt

    #Expand adjacent matrix
    for x in V:
        if x >= len(Adj):
            k_start = len(Adj)
            for k in range(k_start, x+1):
                Adj = np.append(Adj, np.zeros((1, Adj.shape[0])), axis=0)
                Adj = np.append(Adj, np.zeros((Adj.shape[0], 1)), axis=1)
                Adj_opt = np.append(Adj_opt, np.zeros((1, Adj_opt.shape[0])), axis=0)
                Adj_opt = np.append(Adj_opt, np.zeros((Adj_opt.shape[0], 1)), axis=1)
        
    for x, y in itertools.permutations(V,2):#Regard x and y as adjacent
        Adj[x,y] = 1
        for z in AllVertices-{x,y}: #Remove the edges that turns out to be non-adjacent
            if D[x,z]+D[z,y] == D[x,y]:
                Adj[x,y] = 0


    isDeleted = 0

    for x in V:
        print(x,CountNeighbor(x))
        if CountNeighbor(x)==0:
            isDeleted = 1
            V_del.add(x)
        elif CountNeighbor(x)==1:
            isDeleted = 1
            for i in range(len(Adj)):
                if Adj[x,i]==1:
                    u = i
            Adj_opt[x,u] = Adj_opt[u,x] = 1
            V_del.add(x)
        elif CountNeighbor(x)==2:
            isDeleted = 1
            Neighbor = set()
            for i in range(len(Adj)):
                if Adj[x,i]==1:
                    Neighbor.add(i)
            if comp(Neighbor,x)[0] == 0:
                u = Neighbor.pop()
                v = Neighbor.pop()
                Adj_opt[x,u] = Adj_opt[u,x] = 1
                Adj_opt[x,v] = Adj_opt[v,x] = 1
                V_del.add(x)

    for x in V_del:
        for y in V.union(V_del):
            Adj[x,y] = Adj[y,x] = 0

    if isDeleted==0:
        for x,y in itertools.permutations(V.difference(V_del),2):
            if Adj[x,y]==1:
                Adj[x,y] = 0
                Adj_opt[x,y] = 1
        V_del = V_del.union(V)
        return

    
    if len(V.difference(V_del))==1:
        x = V.difference(V_del).pop()
        V_del.add(x)
    
    if len(V.difference(V_del))==2:
        x,y = V.difference(V_del)
        if Adj[x,y]==1:
            Adj[x,y] = 0
            Adj_opt[x,y] = 1
        V_del = V_del.union({x,y})

    if len(V.difference(V_del))>2:
        Compactification(V.difference(V_del))


#Function Slack
#Find all verticies
def Slack(V):
    global D
    global V_del
    global AllVertices
    global Slack_del
    
    V_tmp = copy.deepcopy(V)

    for x in V_tmp:
        if CountOptNeighbor(x)<=2:
            Slack_del.add(x)
        elif CountOptNeighbor(x)>=3:
            Slack_del.add(x)
            Neighbor = set()
            for i in range(len(Adj)):
                if Adj_opt[x,i]==1:
                    Neighbor.add(i)
            V_del = set()
            Neighbor.add(x)
            for x,y in itertools.permutations(Neighbor,2):
                Adj[x,y] = 1
                Adj_opt[x,y] = 0
            Compactification(Neighbor)

    if len(AllVertices.difference(Slack_del))==1:
        x = AllVertices.difference(Slack_del).pop()
        Slack_del.add(x)
    
    if len(AllVertices.difference(Slack_del))==2:
        x,y = AllVertices.difference(Slack_del)
        if Adj[x,y]==1:
            Adj[x,y] = 0
            Adj_opt[x,y] = 1
        Slack_del = Slack_del.union({x,y})

    if len(AllVertices.difference(Slack_del))>2:
        Slack(AllVertices.difference(Slack_del))



#Input
mode = input("csv or stdin: ")
if mode == "csv":#csv mode
    filename = input("File Name: ")
    n, data_list = readcsv(filename)
    D = np.array(data_list)
elif mode == "stdin":#standard input mode
    print("Input n:",end=" ")
    n = int(input()) #Input size of distance matrix
    print("Input Distance Matrix: ")
    D_list = [list(map(float, input().split())) for i in range(n)] #2行目以降標準入力
    D = np.array(D_list)
else:
    print("erorr")
    sys.exit()


X = set(range(n))
VertexNum = n-1

#Construction
Adj = np.array([[0]])
Adj_opt = np.array([[0]])
V_del = set()
AllVertices = copy.deepcopy(X)
Compactification(X)

#Slack
Slack_del = set()
Slack(AllVertices)

#Define the graph
G = nx.Graph(Adj_opt)
for x,y in itertools.permutations(AllVertices,2):
    if Adj_opt[x][y]==1:
        G[x][y]['weight'] = D[x,y]

#Output the graph
edge_labels = nx.get_edge_attributes(G, 'weight')
#If considering weighted plotting, remove the lower three lines
for x,y in itertools.permutations(AllVertices,2):
    if Adj_opt[x][y]==1:
        G[x][y]['weight'] = 1
pos = nx.kamada_kawai_layout(G)
node_color = ['black' if x < n else 'white' for x in AllVertices]
nx.draw_networkx(G, pos=pos, with_labels=False, node_color=node_color, edgecolors='black' , node_size=20)
nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels, font_size=8)
plt.show()