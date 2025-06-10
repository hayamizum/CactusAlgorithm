import numpy as np
import itertools
import os
import csv
import networkx as nx
import matplotlib.pyplot as plt
from typing import Set, Tuple, List, Dict, Optional

class CactusGraphClaude:
    """
    Graph construction algorithm for finite metric spaces
    
    This implementation follows the two-phase algorithm:
    1. Full Compactification: Add auxiliary vertices for slack vertices
    2. Topological Pruning: Identify vertices with degree >= 3 for next iteration
    """
    
    def __init__(self, distance_matrix: np.ndarray):
        """
        Initialize with distance matrix
        
        Args:
            distance_matrix: n×n distance matrix
        """
        self.initial_matrix = distance_matrix.copy()
        self.n = len(distance_matrix)
        self.initial_vertices = set(range(self.n))
        self.max_vertex_id = self.n - 1  # Track maximum vertex ID ever seen
    
    def _matrix_to_dict(self, V: Set[int], matrix: np.ndarray, vertex_to_index: Dict[int, int]) -> Dict[int, Dict[int, float]]:
        """Convert matrix to dict of dicts format"""
        D = {}
        for v in V:
            D[v] = {}
            for u in V:
                v_idx = vertex_to_index[v]
                u_idx = vertex_to_index[u]
                D[v][u] = matrix[v_idx, u_idx]
        return D
    
    def _dict_to_matrix(self, V: Set[int], D: Dict[int, Dict[int, float]]) -> Tuple[np.ndarray, Dict[int, int]]:
        """Convert dict of dicts to matrix format"""
        sorted_vertices = sorted(V)
        vertex_to_index = {v: i for i, v in enumerate(sorted_vertices)}
        n = len(sorted_vertices)
        matrix = np.zeros((n, n))
        
        for i, v in enumerate(sorted_vertices):
            for j, u in enumerate(sorted_vertices):
                matrix[i, j] = D[v][u]
        
        return matrix, vertex_to_index
    
    def _compute_compactification_index(self, v: int, V: Set[int], D: Dict[int, Dict[int, float]]) -> Tuple[float, Optional[Tuple[int, int]]]:
        """Compute compactification index α(v) and the realizing pair"""
        if len(V) < 3:
            return 0.0, None
            
        other_vertices = V - {v}
        if len(other_vertices) < 2:
            return 0.0, None
            
        min_val = float('inf')
        min_pair = None
        
        for u, w in itertools.combinations(other_vertices, 2):
            val = (D[v][u] + D[v][w] - D[u][w]) / 2
            if val < min_val:
                min_val = val
                min_pair = (u, w)
        
        alpha = max(0, min_val)
        return alpha, min_pair if alpha > 1e-9 else None
    
    def _full_compactification(self, V_in: Set[int], D_in: Dict[int, Dict[int, float]]) -> Tuple[Set[int], Dict[int, Dict[int, float]]]:
        """
        Full Compactification: Process each vertex in V_in exactly once
        """
        print(f"\n=== Full Compactification ===")
        print(f"Input vertices: {sorted(V_in)}")
        
        V_out = V_in.copy()
        D_out = {}
        for v in V_out:
            D_out[v] = D_in[v].copy()
        
        # Process each vertex in V_in exactly once
        for v in V_in:
            # Skip if vertex was merged away
            if v not in V_out:
                print(f"  Vertex {v} was merged away, skipping")
                continue
            
            # Compute compactification index
            alpha_v, pair = self._compute_compactification_index(v, V_out, D_out)
            
            if pair is None or alpha_v <= 1e-9:
                continue
            
            u, w = pair
            self.max_vertex_id += 1
            z = self.max_vertex_id
            
            print(f"SLACK VERTEX FOUND: v={v}, alpha_v={alpha_v:.6f}, pair=({u},{w})")
            print(f"  Adding auxiliary vertex z={z}")
            
            # Add new vertex z to V_out and D_out
            V_out.add(z)
            D_out[z] = {}
            
            # Compute distances following cactus_old.py logic
            # Distance from original vertex v to auxiliary vertex z is alpha_v
            D_out[v][z] = D_out[z][v] = alpha_v
            
            # Distances from u and w to z
            D_out[u][z] = D_out[z][u] = max(0, D_out[u][v] - alpha_v)
            D_out[w][z] = D_out[z][w] = max(0, D_out[w][v] - alpha_v)
            
            print(f"  D_out[{u}][{z}] = {D_out[u][z]:.6f}")
            print(f"  D_out[{w}][{z}] = {D_out[w][z]:.6f}")
            
            # Distance from z to itself
            D_out[z][z] = 0.0
            
            # For other vertices x, use max formula from cactus_old.py
            for x in V_out - {v, u, w, z}:
                dist_xz = max(D_out[x][t] - D_out[t][z] for t in {v, u, w})
                dist_xz = max(0, dist_xz)  # Ensure non-negative
                D_out[x][z] = D_out[z][x] = dist_xz
                print(f"  D_out[{x}][{z}] = {dist_xz:.6f}")
            
            # Vertex identification: only check vertices involved in compactification
            # Following cactus_old.py logic - check {v, u, w} for distance 0 to z
            vertices_to_check = {v, u, w}
            merged = False
            
            for a in vertices_to_check:
                if abs(D_out[a][z]) < 1e-9:
                    print(f"  MERGING: auxiliary vertex {z} into {a} (distance = {D_out[a][z]:.6f})")
                    
                    # Remove z from V_out and D_out
                    V_out.remove(z)
                    del D_out[z]
                    for x in V_out:
                        if z in D_out[x]:
                            del D_out[x][z]
                    merged = True
                    break
            
            # If z wasn't merged, check other vertices for distance 0 to z
            if not merged:
                for a in V_out - vertices_to_check - {z}:
                    if abs(D_out[a][z]) < 1e-9:
                        print(f"  MERGING: auxiliary vertex {z} into {a} (distance = {D_out[a][z]:.6f})")
                        
                        # Remove z from V_out and D_out
                        V_out.remove(z)
                        del D_out[z]
                        for x in V_out:
                            if z in D_out[x]:
                                del D_out[x][z]
                        break
            
            if z in V_out:
                print(f"  Auxiliary vertex {z} preserved in V_out")
        
        print(f"Output vertices: {sorted(V_out)}")
        return V_out, D_out
    
    def _compute_non_redundant_edges(self, V: Set[int], D: Dict[int, Dict[int, float]]) -> Set[Tuple[int, int]]:
        """Compute non-redundant edges"""
        edges = set()
        
        for i, j in itertools.combinations(V, 2):
            is_redundant = False
            for k in V - {i, j}:
                if D[i][k] + D[k][j] <= D[i][j] + 1e-9:  # Use epsilon tolerance
                    is_redundant = True
                    break
            
            if not is_redundant:
                edges.add(tuple(sorted((i, j))))
        
        return edges
    
    def _topological_pruning(self, V_in: Set[int], D_in: Dict[int, Dict[int, float]]) -> Set[int]:
        """
        Topological Pruning: Return vertices with degree >= 3 in non-redundant edge graph
        This determines which vertices participate in the next iteration
        """
        print(f"\n=== Topological Pruning ===")
        print(f"Input vertices: {sorted(V_in)}")
        
        # Compute non-redundant edges
        edges = self._compute_non_redundant_edges(V_in, D_in)
        
        if not edges:
            print("No non-redundant edges found")
            return set()
        
        print(f"Non-redundant edges: {len(edges)}")
        
        # Compute vertex degrees
        degrees = {v: 0 for v in V_in}
        for u, v in edges:
            degrees[u] += 1
            degrees[v] += 1
        
        print("Vertex degrees:")
        for v in sorted(V_in):
            print(f"  vertex {v}: degree {degrees[v]}")
        
        # Keep vertices with degree >= 3
        V_out = {v for v in V_in if degrees[v] >= 3}
        
        print(f"Kept vertices (degree >= 3): {sorted(V_out)}")
        print(f"Removed vertices (degree < 3): {sorted(V_in - V_out)}")
        
        return V_out
    
    def _build_final_graph(self, V_final: Set[int], D_final: Dict[int, Dict[int, float]]) -> List[Tuple[int, int, float]]:
        """
        Build final graph from complete graph minus redundant edges
        """
        print(f"\n=== Building Final Graph ===")
        print(f"Final vertices: {sorted(V_final)}")
        
        edges = []
        non_redundant_edges = self._compute_non_redundant_edges(V_final, D_final)
        
        for i, j in non_redundant_edges:
            weight = D_final[i][j]
            edges.append((i, j, weight))
        
        print(f"Final edges: {len(edges)}")
        for u, v, w in edges:
            print(f"  ({u},{v}) weight={w:.6f}")
        
        return edges
    
    def compute(self) -> Tuple[List[Tuple[int, int, float]], Set[int], Dict[int, Dict[int, float]]]:
        """
        Main algorithm implementation
        
        Returns:
            Tuple containing:
            - Edge list: List of (u, v, weight) tuples
            - Final vertex set
            - Final distance matrix as dict of dicts
        """
        # Initialize
        V_current = self.initial_vertices.copy()
        vertex_to_index = {v: v for v in V_current}
        D_current = self._matrix_to_dict(V_current, self.initial_matrix, vertex_to_index)
        
        # Store all vertices and distances for final graph construction
        final_V = V_current.copy()
        final_D = D_current.copy()
        
        iteration = 0
        while True:
            print(f"\n{'='*50}")
            print(f"ITERATION {iteration}")
            print(f"{'='*50}")
            
            V_previous = V_current.copy()
            
            # Phase 1: Full Compactification
            V_compacted, D_compacted = self._full_compactification(V_current, D_current)
            
            # Update final graph construction data only if new vertices were added
            if len(V_compacted) > len(final_V):
                final_V = V_compacted.copy()
                final_D = D_compacted.copy()
                print(f"Updated final graph data: {len(V_compacted)} vertices")
            
            # Phase 2: Topological Pruning
            V_current = self._topological_pruning(V_compacted, D_compacted)
            
            # Update D_current for next iteration (only keep distances between V_current vertices)
            if V_current:
                D_current = {}
                for v in V_current:
                    D_current[v] = {}
                    for u in V_current:
                        D_current[v][u] = D_compacted[v][u]
            
            print(f"\nIteration {iteration} summary:")
            print(f"  V_previous: {sorted(V_previous)}")
            print(f"  V_compacted: {sorted(V_compacted)}")
            print(f"  V_current: {sorted(V_current)}")
            
            # Termination condition: V_current == V_previous
            if V_current == V_previous:
                print(f"\nTerminating: V_current == V_previous")
                break
            
            iteration += 1
        
        # Build final graph using vertices from the iteration with maximum vertices
        edges = self._build_final_graph(final_V, final_D)
        
        print(f"\n=== FINAL RESULT ===")
        print(f"Total vertices: {len(final_V)}")
        print(f"Total edges: {len(edges)}")
        print(f"Initial vertices: {sorted(self.initial_vertices)}")
        print(f"Final vertices: {sorted(final_V)}")
        print(f"Auxiliary vertices: {sorted(final_V - self.initial_vertices)}")
        
        return edges, final_V, final_D

def read_distance_matrix(filename: str) -> Tuple[int, np.ndarray]:
    """Read distance matrix from CSV file"""
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
    """Save graph to file"""
    base_name = os.path.splitext(input_filename)[0]
    output_filename = f"{base_name}_output_claude.txt"
    
    # Sort edges (u < v)
    sorted_edges = sorted([tuple(sorted(e[:2])) + e[2:] for e in graph_edges])
    
    with open(output_filename, 'w') as f:
        f.write("# u, v, weight\n")
        for u, v, weight in sorted_edges:
            f.write(f"{u}, {v}, {weight}\n")
    
    print(f"Graph saved to {output_filename}")

if __name__ == '__main__':
    # User input processing
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
    
    # Compute optimal graph
    print("Processing...")
    cactus = CactusGraphClaude(D)
    edges, vertices, final_D = cactus.compute()
    
    # Save results
    save_graph_to_file(edges, filename)
    print("Done.")
    
    # Graph visualization
    print("Visualizing graph...")
    G = nx.Graph()
    for u, v, weight in edges:
        G.add_edge(u, v, weight=weight)

    # Save actual weights as labels
    edge_labels = nx.get_edge_attributes(G, 'weight')
    
    # Faithfully reproduce cactus_old.py logic:
    # Override all edge weights to 1 before layout calculation
    for u, v in G.edges():
        G[u][v]['weight'] = 1
    
    # Layout is calculated with weight=1 graph
    pos = nx.kamada_kawai_layout(G)
    
    # Color vertices (original vertices black, auxiliary vertices white)
    node_color = ['black' if x < n else 'white' for x in G.nodes()]
    
    nx.draw_networkx(G, pos=pos, with_labels=False, node_color=node_color, 
                    edgecolors='black', node_size=20)
    # Draw edge labels with actual saved weights
    nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels, font_size=8)
    plt.show()