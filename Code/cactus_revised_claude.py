import numpy as np
import itertools
import os
import csv
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
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
        self.vertex_num = self.n - 1  # Track current vertex number (like VertexNum in cactus_old.py)
        self.computation_log = []  # Store computation process log
    
    def _log(self, message: str):
        """Print message and store in log"""
        print(message)
        self.computation_log.append(message)
    
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
        self._log(f"\n=== Full Compactification ===")
        self._log(f"Input vertices: {sorted(V_in)}")
        
        V_out = V_in.copy()
        D_out = {}
        for v in V_out:
            D_out[v] = D_in[v].copy()
        
        # Track auxiliary vertices added in this iteration
        V_aux = set()
        
        # Process each vertex in V_in exactly once
        for v in V_in:
            # Skip if vertex was merged away
            if v not in V_out:
                self._log(f"  Vertex {v} was merged away, skipping")
                continue
            
            # Compute compactification index
            alpha_v, pair = self._compute_compactification_index(v, V_out, D_out)
            
            if pair is None or alpha_v <= 1e-9:
                continue
            
            u, w = pair
            self.vertex_num += 1
            z = self.vertex_num
            
            self._log(f"SLACK VERTEX FOUND: v={v}, alpha_v={alpha_v:.6f}, pair=({u},{w})")
            self._log(f"  Adding auxiliary vertex z={z}")
            
            # Add new vertex z to V_out and D_out
            V_out.add(z)
            D_out[z] = {}
            
            # Compute distances following cactus_old.py logic
            # Distance from original vertex v to auxiliary vertex z is alpha_v
            D_out[v][z] = D_out[z][v] = round(alpha_v, 5)
            
            # Distances from u and w to z
            D_out[u][z] = D_out[z][u] = max(0, round(D_out[u][v] - alpha_v, 5))
            D_out[w][z] = D_out[z][w] = max(0, round(D_out[w][v] - alpha_v, 5))
            
            self._log(f"  D_out[{u}][{z}] = {D_out[u][z]:.6f}")
            self._log(f"  D_out[{w}][{z}] = {D_out[w][z]:.6f}")
            
            # Distance from z to itself
            D_out[z][z] = 0.0
            
            # For other vertices a, use max formula from cactus_old.py
            # D[a,aux] = max(D[a,t] - D[t,aux] for t in {x,y,z})
            # where D[t,aux] has already been computed above
            for a in V_out - {v, u, w, z}:
                candidates = []
                # t = v (x in cactus_old.py): D[a,v] - D[v,z] = D[a,v] - alpha_v
                candidates.append(D_out[a][v] - D_out[v][z])
                # t = u (y in cactus_old.py): D[a,u] - D[u,z]
                candidates.append(D_out[a][u] - D_out[u][z])
                # t = w (z in cactus_old.py): D[a,w] - D[w,z]
                candidates.append(D_out[a][w] - D_out[w][z])
                
                dist_az = max(0, max(candidates))
                D_out[a][z] = D_out[z][a] = dist_az
                self._log(f"  D_out[{a}][{z}] = {dist_az:.6f}")
            
            # Vertex identification: only check vertices involved in compactification
            # Following cactus_old.py logic - check {v, u, w} for distance 0 to z
            vertices_to_check = {v, u, w}
            merged = False
            
            for a in vertices_to_check:
                if abs(D_out[a][z]) < 1e-9:
                    self._log(f"  MERGING: auxiliary vertex {z} into {a} (distance = {D_out[a][z]:.6f})")
                    
                    # Remove z from V_out and D_out
                    V_out.remove(z)
                    del D_out[z]
                    for x in V_out:
                        if z in D_out[x]:
                            del D_out[x][z]
                    self.vertex_num -= 1  # Decrement vertex counter like cactus_old.py
                    merged = True
                    break
            
            # If z wasn't merged, check other vertices for distance 0 to z
            # Following cactus_old.py: only check vertices in V_in.union(V_aux)
            if not merged:
                check_vertices = (V_in.union(V_aux)) & V_out - vertices_to_check - {z}
                for a in check_vertices:
                    if abs(D_out[a][z]) < 1e-9:
                        self._log(f"  MERGING: auxiliary vertex {z} into {a} (distance = {D_out[a][z]:.6f})")
                        
                        # Remove z from V_out and D_out
                        V_out.remove(z)
                        del D_out[z]
                        for x in V_out:
                            if z in D_out[x]:
                                del D_out[x][z]
                        self.vertex_num -= 1  # Decrement vertex counter like cactus_old.py
                        merged = True
                        break
            
            if z in V_out:
                self._log(f"  Auxiliary vertex {z} preserved in V_out")
                V_aux.add(z)  # Track this auxiliary vertex
        
        self._log(f"Output vertices: {sorted(V_out)}")
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
        self._log(f"\n=== Topological Pruning ===")
        self._log(f"Input vertices: {sorted(V_in)}")
        
        # Compute non-redundant edges
        edges = self._compute_non_redundant_edges(V_in, D_in)
        
        if not edges:
            self._log("No non-redundant edges found")
            return set()
        
        self._log(f"Non-redundant edges: {len(edges)}")
        
        # Compute vertex degrees
        degrees = {v: 0 for v in V_in}
        for u, v in edges:
            degrees[u] += 1
            degrees[v] += 1
        
        self._log("Vertex degrees:")
        for v in sorted(V_in):
            self._log(f"  vertex {v}: degree {degrees[v]}")
        
        # Keep vertices with degree >= 3
        V_out = {v for v in V_in if degrees[v] >= 3}
        
        self._log(f"Kept vertices (degree >= 3): {sorted(V_out)}")
        self._log(f"Removed vertices (degree < 3): {sorted(V_in - V_out)}")
        
        return V_out
    
    def _build_final_graph(self, V_final: Set[int], D_final: Dict[int, Dict[int, float]]) -> List[Tuple[int, int, float]]:
        """
        Build final graph from complete graph minus redundant edges
        """
        self._log(f"\n=== Building Final Graph ===")
        self._log(f"Final vertices: {sorted(V_final)}")
        
        edges = []
        non_redundant_edges = self._compute_non_redundant_edges(V_final, D_final)
        
        for i, j in non_redundant_edges:
            weight = D_final[i][j]
            edges.append((i, j, weight))
        
        self._log(f"Final edges: {len(edges)}")
        for u, v, w in edges:
            self._log(f"  ({u},{v}) weight={w:.6f}")
        
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
            self._log(f"\n{'='*50}")
            self._log(f"ITERATION {iteration}")
            self._log(f"{'='*50}")
            
            V_previous = V_current.copy()
            
            # Phase 1: Full Compactification
            V_compacted, D_compacted = self._full_compactification(V_current, D_current)
            
            # Update final graph construction data only if new vertices were added
            if len(V_compacted) > len(final_V):
                final_V = V_compacted.copy()
                final_D = D_compacted.copy()
                self._log(f"Updated final graph data: {len(V_compacted)} vertices")
            
            # Phase 2: Topological Pruning
            V_current = self._topological_pruning(V_compacted, D_compacted)
            
            # Update D_current for next iteration (only keep distances between V_current vertices)
            if V_current:
                D_current = {}
                for v in V_current:
                    D_current[v] = {}
                    for u in V_current:
                        D_current[v][u] = D_compacted[v][u]
            
            self._log(f"\nIteration {iteration} summary:")
            self._log(f"  V_previous: {sorted(V_previous)}")
            self._log(f"  V_compacted: {sorted(V_compacted)}")
            self._log(f"  V_current: {sorted(V_current)}")
            
            # Termination condition: V_current == V_previous
            if V_current == V_previous:
                self._log(f"\nTerminating: V_current == V_previous")
                break
            
            iteration += 1
        
        # Build final graph using vertices from the iteration with maximum vertices
        edges = self._build_final_graph(final_V, final_D)
        
        self._log(f"\n=== FINAL RESULT ===")
        self._log(f"Total vertices: {len(final_V)}")
        self._log(f"Total edges: {len(edges)}")
        self._log(f"Initial vertices: {sorted(self.initial_vertices)}")
        self._log(f"Final vertices: {sorted(final_V)}")
        self._log(f"Auxiliary vertices: {sorted(final_V - self.initial_vertices)}")
        
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

def save_graph_to_file(graph_edges: List[Tuple[int, int, float]], input_filename: str, computation_log: List[str]) -> None:
    """Save graph and computation log to file"""
    base_name = os.path.splitext(input_filename)[0]
    output_filename = f"{base_name}_output_claude.txt"
    
    # Sort edges (u < v)
    sorted_edges = sorted([tuple(sorted(e[:2])) + e[2:] for e in graph_edges])
    
    with open(output_filename, 'w') as f:
        # Write edge list first
        f.write("# u, v, weight\n")
        for u, v, weight in sorted_edges:
            f.write(f"{u}, {v}, {weight}\n")
        
        # Write computation process log
        f.write("\n# COMPUTATION PROCESS LOG\n")
        f.write("# " + "="*48 + "\n")
        for line in computation_log:
            # Handle multi-line entries by splitting on newlines
            for subline in line.split('\n'):
                if subline:  # Skip empty lines from split
                    f.write(f"# {subline}\n")
    
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
    save_graph_to_file(edges, filename, cactus.computation_log)
    print("Done.")
    
    # Verification: Check if computed graph realizes input distance matrix
    print("\nVerifying graph realizes input distance matrix...")
    G_verify = nx.Graph()
    for u, v, weight in edges:
        G_verify.add_edge(u, v, weight=weight)
    
    # Create labels for all vertices in order
    labels = list(range(len(G_verify)))
    # Compute all-pairs shortest distances in the constructed graph
    distances = dict(nx.all_pairs_dijkstra_path_length(G_verify, weight='weight'))
    # Convert to distance matrix (pandas DataFrame)
    distance_matrix = pd.DataFrame(distances).transpose().fillna(float('inf'))
    # Reorder distance matrix by label order
    distance_matrix = distance_matrix[labels].loc[labels]
    # Convert to numpy array
    distance_matrix_np = distance_matrix.values
    # Trim to original n×n matrix (compare only original vertices)
    distance_matrix_trimmed = distance_matrix_np[:n, :n]
    D_original = D[:n, :n]
    # Compare matrices (with rounding to handle floating point precision)
    computed_rounded = np.round(distance_matrix_trimmed, 5)
    original_rounded = np.round(D_original, 5)
    verification_result = np.array_equal(computed_rounded, original_rounded)
    
    if verification_result:
        print("True")
        cactus._log("=== VERIFICATION: True ===")
    else:
        print("False")
        cactus._log("=== VERIFICATION: False ===")
        
        # Find and report mismatched distances
        print("\nMismatched distances:")
        cactus._log("Mismatched distances:")
        mismatch_count = 0
        for i in range(n):
            for j in range(i+1, n):  # Only check upper triangle to avoid duplicates
                original_dist = original_rounded[i, j]
                computed_dist = computed_rounded[i, j]
                if original_dist != computed_dist:
                    mismatch_count += 1
                    msg = f"  Vertices ({i},{j}): Expected {original_dist}, Got {computed_dist}"
                    print(msg)
                    cactus._log(msg)
        
        total_msg = f"Total mismatched pairs: {mismatch_count} out of {n*(n-1)//2} pairs"
        print(total_msg)
        cactus._log(total_msg)
    
    # Update the saved file with verification result
    save_graph_to_file(edges, filename, cactus.computation_log)
    
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
    
    # Draw graph with vertex labels
    nx.draw_networkx(G, pos=pos, with_labels=True, node_color=node_color, 
                    edgecolors='black', node_size=300, font_size=10, font_color='red')
    # Draw edge labels with actual saved weights
    nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels, font_size=8)
    plt.title("Final Graph with Vertex Labels\n(Black: Original vertices, White: Auxiliary vertices)")
    plt.show()