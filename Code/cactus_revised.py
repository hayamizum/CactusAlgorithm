import numpy as np
import itertools
import pandas as pd
import os
import csv
import networkx as nx
import matplotlib.pyplot as plt

# ==================================================================
# データ読み込み & 保存
# ==================================================================

def read_distance_matrix(filename):
    """
    CSVファイルから距離行列を読み込む。
    cactus_old.py の readcsv と同等の機能。
    """
    if not filename.endswith('.csv'):
        filename += '.csv'
        
    with open(filename, 'r', newline='', encoding='utf-8-sig') as file:
        reader = csv.reader(file)
        n = int(next(reader)[0])
        
        distance_matrix_list = []
        for row in reader:
            distance_matrix_list.append([float(value) for value in row])

    V = set(range(n))
    D = {i: {j: distance_matrix_list[i][j] for j in range(n)} for i in range(n)}
    
    return V, D

def save_graph_to_file(graph_edges, input_filename):
    """
    グラフをソート済みの辺リストとしてファイルに保存する。
    cactus_old.py の save_graph_to_file と同等の機能。
    出力ファイル名は input_filename に基づいて生成される。
    """
    base_name = os.path.splitext(input_filename)[0]
    output_filename = f"{base_name}_output_revised.txt"

    # 辺リストを u, v の順でソート (u < v)
    sorted_edges = sorted([tuple(sorted(e[:2])) + e[2:] for e in graph_edges])

    with open(output_filename, 'w') as f:
        f.write("# u, v, weight\n")
        for u, v, weight in sorted_edges:
            f.write(f"{u}, {v}, {weight}\n")
    
    print(f"Graph saved to {output_filename}")

# ==================================================================
# アルゴリズム本体
# ==================================================================

def full_compactification(V_in, D_in):
    """
    疑似コードの Full_Compactification 関数に対応。
    """
    V_out = V_in.copy()
    D_out = {k: v.copy() for k, v in D_in.items()}

    # 新しい頂点を識別するための一意なカウンター
    # 既存の頂点が整数であると仮定し、その最大値から開始する
    # もし頂点名が文字列の場合は、別の方法が必要
    new_vertex_counter = max(list(V_out) + [-1]) + 1

    # ループ中に V_out が変化しても、ループの対象は元の V_in のみ
    for v in list(V_in):
        if v not in V_out: # ループ中に同一視されて消えた場合はスキップ
            continue

        # --- compactification index α(v) の計算 ---
        alpha_v = float('inf')
        compactification_pair = None
        
        other_vertices = V_out - {v}
        if len(other_vertices) < 2:
            continue

        for u, w in itertools.combinations(other_vertices, 2):
            val = (D_out[v][u] + D_out[v][w] - D_out[u][w]) / 2
            if val < alpha_v:
                alpha_v = val
                compactification_pair = (u, w)
        
        # --- 補助点の追加 ---
        # cactus_old.pyでは四捨五入が多用されているが、ここではまず厳密に計算する
        if alpha_v > 1e-9: # 浮動小数点数の誤差を考慮
            u, w = compactification_pair
            z = new_vertex_counter
            new_vertex_counter += 1

            V_out.add(z)
            D_out[z] = {}
            
            # 補助点zと既存の点xとの距離を定義
            # d(x,z) = (d(x,u) + d(x,w) - d(u,w)) / 2
            # この計算には、新しく追加されたz自身も含まれる
            for x in list(V_out):
                # xがu,w,zのいずれかの場合の距離を正しく取得する
                dist_xu = D_out.get(x, {}).get(u, D_out.get(u, {}).get(x))
                dist_xw = D_out.get(x, {}).get(w, D_out.get(w, {}).get(x))
                dist_uw = D_out.get(u, {}).get(w, D_out.get(w, {}).get(u))

                dist_xz = (dist_xu + dist_xw - dist_uw) / 2
                D_out[x][z] = dist_xz
                D_out[z][x] = dist_xz

            # --- 頂点同一視 ---
            while True:
                found_zero_pair = False
                pair_to_merge = None
                
                # 距離が0に近いペアを探す
                for i, j in itertools.combinations(V_out, 2):
                    if abs(D_out[i][j]) < 1e-9: # 浮動小数点数の誤差を考慮
                        pair_to_merge = tuple(sorted((i, j)))
                        found_zero_pair = True
                        break
                
                if found_zero_pair:
                    i, j = pair_to_merge # i < j
                    
                    # j を i に統合する (j を削除)
                    V_out.remove(j)
                    
                    # D_out から j の行と列を削除
                    del D_out[j]
                    for k in D_out:
                        del D_out[k][j]
                else:
                    break # 同一視するペアがなくなったら終了
    
    return V_out, D_out

def topological_pruning(V_in, D_in, V_initial):
    """
    【この関数は不要になりました】
    疑似コードの Topological_Pruning 関数に対応。
    V_initial に含まれる頂点は削除しないように変更。
    """
    # --- 冗長でない辺集合 E_non_redundant の構築 ---
    E_non_redundant = set()
    for i, j in itertools.combinations(V_in, 2):
        is_redundant = False
        # 辺(i,j)が他の頂点kを経由して冗長かどうかを判定
        for k in V_in - {i, j}:
            # d(i,k) + d(k,j) <= d(i,j)
            # 浮動小数点数の誤差を考慮
            if D_in[i][k] + D_in[k][j] <= D_in[i][j] + 1e-9:
                is_redundant = True
                break
        
        if not is_redundant:
            E_non_redundant.add(tuple(sorted((i, j))))

    # --- 次数が3以上の頂点を抽出 ---
    V_out = set()
    if not E_non_redundant:
        return V_out

    # 各頂点の次数を計算
    degrees = {v: 0 for v in V_in}
    for u, v in E_non_redundant:
        degrees[u] += 1
        degrees[v] += 1
    
    for v, degree in degrees.items():
        # 次数が3以上、または元の頂点であれば維持する
        if degree >= 3 or v in V_initial:
            V_out.add(v)
            
    return V_out

def build_final_graph(V_final, D_final):
    """
    疑似コードの Build_Final_Graph 関数に対応。
    """
    final_edges = []
    for i, j in itertools.combinations(V_final, 2):
        is_redundant = False
        for k in V_final - {i, j}:
            # d(i,k) + d(k,j) <= d(i,j)
            if D_final[i][k] + D_final[k][j] <= D_final[i][j] + 1e-9:
                is_redundant = True
                break
        
        if not is_redundant:
            final_edges.append((i, j, D_final[i][j]))
            
    return final_edges

def main(V_initial, D_initial):
    """
    疑似コードの Main 関数に対応。
    中間的なPruningは不要であるという発見に基づき、ロジックを修正。
    """
    V_current = V_initial.copy()
    D_current = D_initial.copy()

    while True:
        V_previous_size = len(V_current)

        # フェーズ1: Full Compactification & 頂点同一視
        # これを、頂点数が増えなくなるまで繰り返す
        V_current, D_current = full_compactification(V_current, D_current)

        # --- 終了条件 ---
        # Compactificationで頂点数に変化がなければループ終了
        if len(V_current) == V_previous_size:
            break
            
    # --- 仕上げ: 最終的なグラフの構築 ---
    # ループで完成した最終的な頂点集合と距離行列を使う
    final_graph_edges = build_final_graph(V_current, D_current)
    
    return final_graph_edges, V_current, D_current

# ==================================================================
# スクリプト実行のエントリポイント
# ==================================================================

if __name__ == '__main__':
    # ユーザーからの入力を受け付ける部分
    # (cactus_old.py と同様に csv / stdin を選択できるようにする)
    mode = input("csv or stdin: ")
    filename = ""
    n_initial = 0

    if mode.lower() in ["csv", "c"]:
        filename = input("File Name: ")
        V, D = read_distance_matrix(filename)
        n_initial = len(V)
    elif mode.lower() in ["stdin", "s"]:
        print("Input n:", end=" ")
        n_initial = int(input())
        print("Input Distance Matrix: ")
        D_list = [list(map(float, input().split())) for i in range(n_initial)]
        V = set(range(n_initial))
        D = {i: {j: D_list[i][j] for j in range(n_initial)} for i in range(n_initial)}
        filename = "stdin_input"
    else:
        print("error: invalid mode")
        exit()
    
    print("Processing...")
    final_graph_edges, V_final, D_final = main(V, D)
    save_graph_to_file(final_graph_edges, filename)
    print("Done.")

    # --- グラフ描画処理 ---
    print("Visualizing graph...")
    G = nx.Graph()
    for u, v, weight in final_graph_edges:
        G.add_edge(u, v, weight=weight)

    # 実際の重みをラベルとして保存
    edge_labels = nx.get_edge_attributes(G, 'weight')
    
    # cactus_old.pyのロジックを忠実に再現:
    # レイアウト計算の前に、全ての辺の重みを1に上書きする
    for u, v in G.edges():
        G[u][v]['weight'] = 1
    
    # これで、重み1のグラフとしてレイアウトが計算される
    pos = nx.kamada_kawai_layout(G)
    
    # 頂点の色分け (元頂点は黒、補助頂点は白)
    node_color = ['black' if x < n_initial else 'white' for x in G.nodes()]
    
    nx.draw_networkx(G, pos=pos, with_labels=False, node_color=node_color, edgecolors='black' , node_size=20)
    # ラベルは保存しておいた実際の重みで描画する
    nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels, font_size=8)
    plt.show()
