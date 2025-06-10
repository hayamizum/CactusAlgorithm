import numpy as np
import itertools
import pandas as pd
import os
import csv
import networkx as nx
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Set, Tuple, List, Optional
import logging

# ロガーの設定
logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter('%(message)s'))
logger.addHandler(handler)

@dataclass
class CompactificationResult:
    """Compactification操作の結果を保持するデータクラス"""
    vertex_set: Set[int]  # 頂点集合
    distance_matrix: np.ndarray  # 距離行列
    auxiliary_vertices: Set[int]  # 追加された補助点の集合
    vertex_to_index: dict  # 頂点IDから行列インデックスへのマッピング
    index_to_vertex: dict  # 行列インデックスから頂点IDへのマッピング

class CactusGraph:
    """
    距離空間からカクタスグラフを構築するクラス
    
    Attributes:
        D (np.ndarray): 距離行列
        n (int): 入力の頂点数
        initial_vertices (Set[int]): 初期頂点集合
        debug (bool): デバッグモードフラグ
    """
    def __init__(self, distance_matrix: np.ndarray, debug: bool = False):
        """
        Args:
            distance_matrix: n×nの距離行列
            debug: デバッグ出力を有効にするかどうか
        """
        self.D = distance_matrix.copy()
        self.n = len(distance_matrix)
        self.initial_vertices = set(range(self.n))
        self.debug = debug
        if debug:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.INFO)
    
    def _log_state(self, phase: str, V: Set[int], D: np.ndarray, vertex_to_index: dict) -> None:
        """現在の状態をログ出力"""
        if not self.debug:
            return
        
        logger.debug(f"\n=== {phase} ===")
        logger.debug(f"Vertices: {sorted(list(V))}")
        logger.debug(f"Vertex to Index mapping: {vertex_to_index}")
        logger.debug("Distance Matrix:")
        with np.printoptions(precision=3, suppress=True):
            logger.debug(str(D))

    def _compute_compactification_index(self, v: int, V: Set[int], D: np.ndarray, vertex_to_index: dict) -> Tuple[float, Optional[Tuple[int, int]]]:
        """頂点vのcompactification indexとその実現ペアを計算"""
        if len(V) < 3:
            return 0.0, None
            
        other_vertices = V - {v}
        min_val = float('inf')
        min_pair = None
        v_idx = vertex_to_index[v]
        
        for u, w in itertools.combinations(other_vertices, 2):
            u_idx = vertex_to_index[u]
            w_idx = vertex_to_index[w]
            val = (D[v_idx,u_idx] + D[v_idx,w_idx] - D[u_idx,w_idx]) / 2
            if val < min_val:
                min_val = val
                min_pair = (u, w)
                
            if self.debug:
                logger.debug(f"  Checking pair ({u},{w}) for vertex {v}: val = {val}")
        
        result = max(0, round(min_val, 5)), min_pair
        if self.debug:
            logger.debug(f"  Compactification index for vertex {v}: {result[0]} with pair {result[1]}")
        return result

    def _full_compactification(self, V: Set[int], D: np.ndarray) -> CompactificationResult:
        """
        Full Compactification操作を実行
        
        Args:
            V: 処理対象の頂点集合
            D: 現在の距離行列
            
        Returns:
            CompactificationResult: 操作結果
        """
        # 頂点IDと行列インデックスのマッピングを初期化
        vertex_to_index = {v: i for i, v in enumerate(sorted(V))}
        index_to_vertex = {i: v for v, i in vertex_to_index.items()}
        
        self._log_state("Starting Full Compactification", V, D, vertex_to_index)
        
        V_out = V.copy()
        D_out = D.copy()
        auxiliary_vertices = set()
        
        # 新しい頂点のIDは現在の最大値+1から開始
        next_vertex_id = max(V) + 1
        
        for v in list(V):  # Vのコピー上でループ
            if v not in V_out:  # 同一視で削除された可能性
                continue
                
            alpha_v, pair = self._compute_compactification_index(v, V_out, D_out, vertex_to_index)
            if alpha_v <= 1e-9 or pair is None:  # 浮動小数点誤差を考慮
                continue
                
            u, w = pair
            z = next_vertex_id
            next_vertex_id += 1
            
            if self.debug:
                logger.debug(f"\nAdding auxiliary vertex {z} for vertex {v}")
                logger.debug(f"  Based on pair ({u},{w}) with alpha = {alpha_v}")
            
            # 距離行列の拡張
            old_size = len(D_out)
            new_D = np.zeros((old_size + 1, old_size + 1))
            new_D[:old_size, :old_size] = D_out
            
            # 新頂点zと他の頂点との距離を計算
            z_idx = old_size  # 新しい頂点の行列インデックス
            v_idx = vertex_to_index[v]
            u_idx = vertex_to_index[u]
            w_idx = vertex_to_index[w]

            # d(v,z), d(u,z), d(w,z) を正しく計算
            dist_vz = alpha_v
            dist_uz = D_out[u_idx, v_idx] - alpha_v
            dist_wz = D_out[w_idx, v_idx] - alpha_v

            new_D[v_idx, z_idx] = new_D[z_idx, v_idx] = dist_vz
            new_D[u_idx, z_idx] = new_D[z_idx, u_idx] = dist_uz
            new_D[w_idx, z_idx] = new_D[z_idx, w_idx] = dist_wz

            # v,u,w以外の頂点xとzの距離を計算
            other_vertices_in_V = V_out - {v, u, w}
            for x in other_vertices_in_V:
                x_idx = vertex_to_index[x]
                dist_xz = (D_out[x_idx,u_idx] + D_out[x_idx,w_idx] - D_out[u_idx,w_idx]) / 2
                new_D[x_idx,z_idx] = new_D[z_idx,x_idx] = dist_xz
            
            new_D[z_idx,z_idx] = 0
            
            # マッピングを更新
            vertex_to_index[z] = z_idx
            index_to_vertex[z_idx] = z
            
            D_out = new_D
            V_out.add(z)
            auxiliary_vertices.add(z)
            
            self._log_state(f"After adding vertex {z}", V_out, D_out, vertex_to_index)
            
            # 距離0の頂点ペアを同一視
            while True:
                zero_pairs = []
                for i, j in itertools.combinations(V_out, 2):
                    i_idx = vertex_to_index[i]
                    j_idx = vertex_to_index[j]
                    if abs(D_out[i_idx,j_idx]) < 1e-9:
                        zero_pairs.append((i, j))
                if not zero_pairs:
                    break
                    
                i, j = min(zero_pairs)  # 決定的な結果のため最小のペアを選択
                if self.debug:
                    logger.debug(f"\nMerging vertices {j} into {i} (distance = {D_out[vertex_to_index[i],vertex_to_index[j]]})")
                
                # jをiに統合（jを削除）
                j_idx = vertex_to_index[j]
                V_out.remove(j)
                if j in auxiliary_vertices:
                    auxiliary_vertices.remove(j)
                
                # j行j列を削除
                mask = np.ones(len(D_out), dtype=bool)
                mask[j_idx] = False
                D_out = D_out[mask][:,mask]
                
                # マッピングを更新
                del vertex_to_index[j]
                del index_to_vertex[j_idx]
                # j_idxより大きいインデックスを持つ頂点のインデックスを更新
                for v in V_out:
                    if vertex_to_index[v] > j_idx:
                        vertex_to_index[v] -= 1
                index_to_vertex = {i: v for v, i in vertex_to_index.items()}
                
                self._log_state("After merging vertices", V_out, D_out, vertex_to_index)
        
        return CompactificationResult(V_out, D_out, auxiliary_vertices, vertex_to_index, index_to_vertex)

    def _compute_non_redundant_edges(self, V: Set[int], D: np.ndarray, vertex_to_index: dict) -> Set[Tuple[int, int]]:
        """冗長でない辺の集合を計算"""
        edges = set()
        for i, j in itertools.combinations(V, 2):
            is_redundant = False
            i_idx = vertex_to_index[i]
            j_idx = vertex_to_index[j]
            for k in V - {i, j}:
                k_idx = vertex_to_index[k]
                if D[i_idx,k_idx] + D[k_idx,j_idx] <= D[i_idx,j_idx] + 1e-9:  # 浮動小数点誤差を考慮
                    is_redundant = True
                    break
            if not is_redundant:
                edges.add(tuple(sorted((i, j))))
                
            if self.debug and not is_redundant:
                logger.debug(f"  Non-redundant edge found: ({i},{j}) with weight {D[i_idx,j_idx]}")
                
        return edges

    def _topological_pruning(self, V: Set[int], D: np.ndarray, vertex_to_index: dict) -> Set[int]:
        """
        次数3以上の頂点および初期頂点を保持
        
        Args:
            V: 頂点集合
            D: 距離行列
            vertex_to_index: 頂点IDから行列インデックスへのマッピング
            
        Returns:
            Set[int]: 保持する頂点の集合
        """
        self._log_state("Starting Topological Pruning", V, D, vertex_to_index)
        
        edges = self._compute_non_redundant_edges(V, D, vertex_to_index)
        if not edges:
            return set()
            
        # 各頂点の次数を計算
        degrees = {v: 0 for v in V}
        for u, v in edges:
            degrees[u] += 1
            degrees[v] += 1
            
        if self.debug:
            logger.debug("\nVertex degrees:")
            for v in sorted(V):
                logger.debug(f"  vertex {v}: degree {degrees[v]}")
        
        # 次数3以上または初期頂点を保持
        result = {v for v in V if degrees[v] >= 3 or v in self.initial_vertices}
        
        if self.debug:
            logger.debug(f"\nKept vertices: {sorted(list(result))}")
        
        return result

    def compute(self) -> Tuple[List[Tuple[int, int, float]], Set[int], np.ndarray]:
        """
        カクタスグラフを計算
        
        Returns:
            Tuple[List[Tuple[int, int, float]], Set[int], np.ndarray]:
                - エッジリスト（各エッジは(始点,終点,重み)のタプル）
                - 最終的な頂点集合
                - 最終的な距離行列
        """
        V_current = self.initial_vertices.copy()
        D_current = self.D.copy()
        vertex_to_index = {v: v for v in V_current}  # 初期状態では頂点IDとインデックスは同じ
        
        iteration = 0
        while True:
            if self.debug:
                logger.debug(f"\n=== Iteration {iteration} ===")
            
            V_previous = V_current.copy()
            
            # フェーズ1: Full Compactification
            result = self._full_compactification(V_current, D_current)
            
            # フェーズ2: Topological Pruning
            V_current = self._topological_pruning(result.vertex_set, result.distance_matrix, result.vertex_to_index)
            
            # 距離行列を更新
            if V_current == V_previous:
                break
                
            # 保持する頂点の行と列のみを抽出
            indices = [result.vertex_to_index[v] for v in sorted(V_current)]
            D_current = result.distance_matrix[np.ix_(indices, indices)]
            vertex_to_index = {v: i for i, v in enumerate(sorted(V_current))}
            
            iteration += 1
        
        # 最終的なグラフの構築
        edges = []
        for i, j in self._compute_non_redundant_edges(result.vertex_set, result.distance_matrix, result.vertex_to_index):
            i_idx = result.vertex_to_index[i]
            j_idx = result.vertex_to_index[j]
            edges.append((i, j, result.distance_matrix[i_idx,j_idx]))
        
        if self.debug:
            logger.debug("\n=== Final Result ===")
            logger.debug(f"Vertices: {sorted(list(result.vertex_set))}")
            logger.debug(f"Edges: {edges}")
        
        return edges, result.vertex_set, result.distance_matrix

def read_distance_matrix(filename: str) -> Tuple[int, np.ndarray]:
    """CSVファイルから距離行列を読み込む"""
    if not filename.endswith('.csv'):
        filename += '.csv'
        
    with open(filename, 'r', newline='', encoding='utf-8-sig') as file:
        reader = csv.reader(file)
        n = int(next(reader)[0])
        
        distance_matrix = []
        for row in reader:
            distance_matrix.append([float(value) for value in row])
            
    return n, np.array(distance_matrix)

def save_graph_to_file(graph_edges: List[Tuple[int, int, float]], input_filename: str) -> None:
    """グラフをファイルに保存"""
    base_name = os.path.splitext(input_filename)[0]
    output_filename = f"{base_name}_output_revised.txt"
    
    # 辺をソート (u < v)
    sorted_edges = sorted([tuple(sorted(e[:2])) + e[2:] for e in graph_edges])
    
    with open(output_filename, 'w') as f:
        f.write("# u, v, weight\n")
        for u, v, weight in sorted_edges:
            f.write(f"{u}, {v}, {weight}\n")
    
    print(f"Graph saved to {output_filename}")

if __name__ == '__main__':
    # ユーザー入力の処理
    mode = input("csv or stdin: ")
    filename = ""
    
    if mode.lower() in ["csv", "c"]:
        filename = input("File Name: ")
        n, D = read_distance_matrix(filename)
    elif mode.lower() in ["stdin", "s"]:
        print("Input n:", end=" ")
        n = int(input())
        print("Input Distance Matrix: ")
        D_list = [list(map(float, input().split())) for _ in range(n)]
        D = np.array(D_list)
        filename = "stdin_input"
    else:
        print("error: invalid mode")
        exit()
    
    # デバッグモードの設定
    debug_mode = input("Enable debug mode? (y/n): ").lower() == 'y'
    
    # カクタスグラフの計算
    print("Processing...")
    cactus = CactusGraph(D, debug=debug_mode)
    edges, vertices, final_D = cactus.compute()
    
    # 結果の保存
    save_graph_to_file(edges, filename)
    print("Done.")
    
    # --- グラフ描画処理 ---
    print("Visualizing graph...")
    G = nx.Graph()
    for u, v, weight in edges:
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
    node_color = ['black' if x < n else 'white' for x in G.nodes()]
    
    nx.draw_networkx(G, pos=pos, with_labels=False, node_color=node_color, edgecolors='black', node_size=20)
    # ラベルは保存しておいた実際の重みで描画する
    nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels, font_size=8)
    plt.show()
