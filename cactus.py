#Constructing Cactus
#
#written by Keita Watanabe, Waseda university
#2023.1.16
#
#
import sys
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import itertools
import copy
import csv
import pandas as pd

#関数 readcsv
def readcsv(filename):
    distance_matrix = []
    with open(filename, 'r', newline='', encoding='utf-8-sig') as file:
        reader = csv.reader(file)
        for i, row in enumerate(reader):
            if i == 0:
                # 1行目は行列の大きさの情報
                n = int(row[0])
                continue
            # 2行目以降のデータを数値としてリストに追加します
            distance_matrix.append([float(value) for value in row])
    return n, np.array(distance_matrix)

#関数 add_csv .csvを付ける
def add_csv(filename):
    return filename if filename.endswith('.csv') else f"{filename}.csv"

#関数 non_negative
def non_negative(value):
    return max(value, 0)

#深さ優先探索
def DFS(graph, v, visited, component):
    visited[v] = True
    component.append(v)
    for i in range(len(graph)):
        if graph[v][i] == 1 and not visited[i]:
            DFS(graph, i, visited, component)

#連結成分
def connected_components(graph):
    num_vertices = len(graph)
    visited = [False] * num_vertices
    components = []

    for v in range(num_vertices):
        if not visited[v]:
            component = []
            DFS(graph, v, visited, component)
            components.append(component)

    return components


#関数comp Compactification indexをV内で計算する．
def comp(V, x):

    global D
    
    #x以外の要素をy,zに初期値として割り当てる
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

    
    return non_negative(round(min_val,5)) , min_comb

#関数Compactification
def Compactification(V):
    global D
    global VertexNum
    global V_del
    global AllVertices
    global Adj
    global Adj_opt
    V_aux = set()

    Adj = np.zeros((len(Adj),len(Adj)))

    if len(V)<=2:
        mkAdj(V)
        return

    for x in V:
        index, (y, z) = comp(V,x)
        if index == 0:
            continue
        else:
            #追加頂点の情報
            VertexNum += 1
            aux = VertexNum
            #距離行列拡張
            D = np.append(D, np.zeros((1, D.shape[0])), axis=0)
            D = np.append(D, np.zeros((D.shape[0], 1)), axis=1)
            #距離行列書き込み(まずcompに関わったもの)
            D[x,aux] = D[aux,x] = index
            D[y,aux] = D[aux,y] = non_negative(round(D[y,x] - D[x,aux],5))
            D[z,aux] = D[aux,z] = non_negative(round(D[z,x] - D[x,aux],5))
            #同一視
            isEqualed = 0
            for a in {x,y,z}:
                if D[a, aux] == 0:
                    VertexNum -= 1
                    D = np.delete(D, aux, axis=0)
                    D = np.delete(D, aux, axis=1)
                    isEqualed = 1
                    break
            if isEqualed == 1:
                continue#この後はやらなくてよい
            #距離行列書き込み(compに関わってないもの)
            for a in AllVertices.union(V_aux)-{x,y,z,aux}:
                D[a,aux] = D[aux,a] = non_negative(round(max(D[a,t] - D[t,aux] for t in {x,y,z}),5))
                #同一視
                if a in V.union(V_aux):
                    if D[a, aux] == 0:
                        VertexNum -= 1
                        D = np.delete(D, aux, axis=0)
                        D = np.delete(D, aux, axis=1)
                        isEqualed = 1
                        break
            if isEqualed == 0:
                V_aux.add(aux)
                AllVertices.add(aux)

    mkAdj(V.union(V_aux))



#関数CountNeighbor
def CountNeighbor(x, Adj):
    NeighborNum = 0
    for i in range(len(Adj)):
        if Adj[x,i]==1:
            NeighborNum += 1
    return NeighborNum

#関数makeAdjacent
def mkAdj(V):
    global D
    global V_del
    global AllVertices
    global Adj
    global Adj_opt
    global Adj_opt_tmp

    #隣接行列拡張
    for x in V:
        if x >= len(Adj):
            k_start = len(Adj)
            for k in range(k_start, x+1):
                Adj = np.append(Adj, np.zeros((1, Adj.shape[0])), axis=0)
                Adj = np.append(Adj, np.zeros((Adj.shape[0], 1)), axis=1)
                Adj_opt = np.append(Adj_opt, np.zeros((1, Adj_opt.shape[0])), axis=0)
                Adj_opt = np.append(Adj_opt, np.zeros((Adj_opt.shape[0], 1)), axis=1)
                Adj_opt_tmp = np.append(Adj_opt_tmp, np.zeros((1, Adj_opt_tmp.shape[0])), axis=0)
                Adj_opt_tmp = np.append(Adj_opt_tmp, np.zeros((Adj_opt_tmp.shape[0], 1)), axis=1)
        
    for x, y in itertools.permutations(V,2): #完全グラフの隣接関係をいれる
        Adj[x,y] = 1
        for z in V-{x,y}: #隣接していないものを切る
            if D[x,z]+D[z,y] == D[x,y] and D[x,z] != 0 and D[z,y] != 0:
                Adj[x,y] = Adj[y,x] = 0


    isDeleted = 0

    for x in V:
        if CountNeighbor(x, Adj)==0:
            isDeleted = 1
            V_del.add(x)
        elif CountNeighbor(x, Adj)==1:
            isDeleted = 1
            for i in range(len(Adj)):
                if Adj[x,i]==1:
                    u = i
            Adj_opt[x,u] = Adj_opt[u,x] = 1
            V_del.add(x)
        elif CountNeighbor(x, Adj)==2:
            isDeleted = 1
            Neighbor = set()
            for i in range(len(Adj)):
                if Adj[x,i]==1:
                    Neighbor.add(i)
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
                Adj[x,y] = Adj[y,x] = 0
                Adj_opt[x,y] = Adj_opt[y,x] = 1
        V_del = V_del.union(V)
        return


    components = connected_components(Adj)
    
    for connect in components:
        connect = set(connect)
        if len(connect)==1:
            x = connect.pop()
            V_del.add(x)
    
        if len(connect)==2:
            x,y = connect
            if Adj[x,y]==1:
                Adj[x,y] = Adj[y,x] = 0
                Adj_opt[x,y] = Adj_opt[y,x] = 1
            V_del = V_del.union({x,y})

        if len(connect)>2:
            Compactification(connect)


#関数SetDelete
def SetDelete(V, delete):
    V_list = list(V.difference(delete))
    V_list.sort()
    for x in V_list:
        V_list[V_list.index(x)] = V_list.index(x)
    return set(V_list)


#関数Slack
def Slack(V):
    global D
    global V_del
    global AllVertices
    global Slack_del
    global Adj
    global Adj_opt
    global Adj_opt_tmp
    global VertexNum
    global loop
    
    V_tmp = copy.deepcopy(V)
    Adj_opt_tmp = copy.deepcopy(Adj_opt)

    loop = loop + 1
    Set_of_Neighbor = set()

    for x in V_tmp:
        if CountNeighbor(x, Adj_opt_tmp)<=1:
            Slack_del.add(x)
            continue
        elif CountNeighbor(x, Adj_opt_tmp)>=2:
            Neighbor = set()
            for i in range(len(Adj_opt_tmp)):
                if Adj_opt_tmp[x,i]==1:
                    Neighbor.add(i)
            Neighbor.add(x)
            Slack_del.add(x)
            Set_of_Neighbor.add(frozenset(Neighbor))
            for y,z in itertools.permutations(Neighbor,2):
                Adj_opt[y,z] = Adj_opt[z,y] = 0

    for Neighbor in Set_of_Neighbor:
        Compactification(Neighbor)

    Delete_set = set()
    for x,y in itertools.combinations(AllVertices,2):
        if D[x,y]==0:
            Delete_set.add(max(x,y))
            for i in range(len(Adj_opt)):
                if i==min(x,y):
                    Adj_opt[min(x,y),i] = 0
                else:
                    Adj_opt[min(x,y),i] = Adj_opt[x,i] or Adj_opt[y,i]
                    Adj_opt[i,min(x,y)] = Adj_opt[i,x] or Adj_opt[i,y]

    Delete_list = list(Delete_set)

    D = np.delete(D, Delete_list, axis=0)
    D = np.delete(D, Delete_list, axis=1)
    Adj = np.delete(Adj, Delete_list, axis=0)
    Adj = np.delete(Adj, Delete_list, axis=1)
    Adj_opt = np.delete(Adj_opt, Delete_list, axis=0)
    Adj_opt = np.delete(Adj_opt, Delete_list, axis=1)

    VertexNum = VertexNum - len(Delete_list)
    AllVertices = SetDelete(AllVertices, Delete_list)
    V = SetDelete(V, Delete_list)
    V_tmp = SetDelete(V_tmp, Delete_list)
    V_del = SetDelete(V_del, Delete_list)
    Slack_del = SetDelete(Slack_del, Delete_list)
    
    if loop==1:
        Slack(AllVertices)
    elif len(AllVertices.difference(Slack_del))>0:
        Slack(AllVertices.difference(Slack_del))



#Input
mode = input("csv or stdin: ")
if mode == "csv" or mode == "c":#ファイル名を標準入力から入力してデータをリストに格納
    filename = input("File Name: ")
    n, data_list = readcsv(add_csv(filename))
    D = np.array(data_list)
elif mode == "stdin" or mode == "s":#標準入力
    print("Input n:",end=" ")
    n = int(input()) #行数の入力
    print("Input Distance Matrix: ")
    D_list = [list(map(float, input().split())) for i in range(n)] #2行目以降標準入力
    D = np.array(D_list)
else:
    print("erorr")
    sys.exit()


X = set(range(n))
VertexNum = n-1

#構成
Adj = np.ones((n,n))
np.fill_diagonal(Adj, 0)
Adj_opt = np.ones((n,n))
np.fill_diagonal(Adj_opt, 0)
Adj_opt_tmp = np.ones((n,n))
np.fill_diagonal(Adj_opt_tmp, 0)
V_del = set()
AllVertices = copy.deepcopy(X)
Slack_del = set()
loop = 0
Slack(X)


#部分的構築により要らない辺がある場合があるため要らない辺を切る
for x, y in itertools.permutations(AllVertices,2): #完全グラフの隣接関係をいれる
    Adj_opt[x,y] = 1
    for z in AllVertices-{x,y}: #隣接していないものを切る
        if D[x,z]+D[z,y] == D[x,y] and D[x,z] != 0 and D[z,y] != 0:
            Adj_opt[x,y] = Adj_opt[y,x] =0

#Graphの定義
G = nx.Graph(Adj_opt)
for x,y in itertools.permutations(AllVertices,2):
    if Adj_opt[x][y]==1:
        G[x][y]['weight'] = round(D[x,y],5)

edge_labels = nx.get_edge_attributes(G, 'weight')


#元データとあっているか
# ラベルを0からの順序で作成する
labels = list(range(len(G)))
# 頂点間の全てのペアに対して最短距離を計算します．
distances = dict(nx.all_pairs_dijkstra_path_length(G, weight='weight'))
# 距離行列（PandasのDataFrame）に変換します．
distance_matrix = pd.DataFrame(distances).transpose().fillna(float('inf'))
# ラベルの並び順に距離行列を並び替えます．
distance_matrix = distance_matrix[labels].loc[labels]
#npに変換
distance_matrix_np = distance_matrix.values
# numpyのarrayに変換されたdistance_matrixをnxnにトリミング
distance_matrix_trimmed = distance_matrix_np[:n, :n]
D_loc = D[:n, :n]
# distance_matrixとDを比較
if np.array_equal(np.round(distance_matrix_trimmed,5), np.round(D_loc,5)):
    print("True")
else:
    print("False")

#重みを考慮した配置にしたい場合は以下の三行をコメントアウト
for x,y in itertools.permutations(AllVertices,2):
    if Adj_opt[x][y]==1:
        G[x][y]['weight'] = 1
pos = nx.kamada_kawai_layout(G)
node_color = ['black' if x < n else 'white' for x in AllVertices]
nx.draw_networkx(G, pos=pos, with_labels=False, node_color=node_color, edgecolors='black' , node_size=20)
nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels, font_size=8)
plt.show()