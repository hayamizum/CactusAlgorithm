# Project Introduction

This project aims at computing a graph realization of a given distance matrix. A distance matrix contains the distances between all object pairs. For example, if a distance matrix of 10 objects is given, then there should be 10×10 elements. A realization is a graph G = (V,E) whose shortest path distances d_G between all pairs of objects in V are the same as given in the distance matrix. In the case when the input distance matrix can be realized by a cactus graph or a fully labeled graph, the program outputs its unique optimal realization (i.e., a realization that minimizes the sum of edge-weights of all edges of the graph).

This project is written in Python and can be run in a Python environment.

## Environment Setup

First, clone this repository to your local machine and access the main directory using the command below:

```bash
git clone https://github.com/hayamizum/CactusAlgorithm.git
cd CactusAlgorithm
```

To run this project, use:

```bash
python cactus.py
```

and follow the instructions in tutorial.

## Prerequisites

To run this project, you may require the following packages:

- sys
- numpy
- networkx
- matplotlib
- itertools
- copy
- csv
- pandas

## Tutorial

When running the code, there are two modes for input dataset: csv or stdin.

For csv mode, it is required to provide the dataset by csv format file. A path of absolute path or relative path is needed for searching the file.

The csv file contains n+1 lines, in which n represents the number of objects, and should follow the criterias below:

- in the first line, the number of objects n should be given, like 5.
- in the following n lines, the distance matrix should be given, with n numbers every line.

As an example, a file called 'manhattan_distances_20x20.csv' which contains L1 distance (Manhattan Distance) matrix of 20 random points in a 20*20 grid is given. By following the instruction below, a graph containing the exact description of distance matrix is shown. 
```
csv or stdin:csv
File Name:manhattan_distances_20x20.csv
```
The graph should be like this (for the simplicity, I manually chose not to show the weights of edges, but the result should contain. Also, the picture is mirror flipped and rotated, but the structure does not change).

<img width="565" alt="Untitled (3)" src="https://github.com/keita1126/CactusAlgorithm/assets/31284538/6a0229ff-2cd6-4aeb-a2c4-7a1297d80ea5">

The original distribution of points is this (you can use generate_random_integer_points.py to generate the random points yourself):
![Untitled (2)](https://github.com/keita1126/CactusAlgorithm/assets/31284538/acefae24-8ed6-47cb-b2e9-9a5dd6a1eaec)






For stdin mode, it is required to input the dataset manually. The first parameter is about the number of objects n. Then the input distance matrix including n×n elements is required, and the data should use white space and Enter to distinguish elements.

For example:

```
Input n: 8
Input Distance Matrix: 
0 1 2 1 1 2 3 2
1 0 1 2 2 1 2 3
2 1 0 1 3 2 1 2
1 2 1 0 2 3 2 1
1 2 3 2 0 1 2 1
2 1 2 3 1 0 1 2
3 2 1 2 2 1 0 1
2 3 2 1 1 2 1 0
```

The result should be a picture of a cube with every edge in weight 1. Also, a result about whether the generated graph is the realization of the given distance matrix is shown as True or False.
![image](https://github.com/keita1126/CactusAlgorithm/assets/31284538/7ed85f02-f9ba-46c3-b171-0a8f417fea2a)

## Special Examples

### For HIV data

The original HIV sequence fasta file has been provided as `hiv-db_gap_strip.fasta`, the `hiv_seq_hamming_distance.csv` is the hamming distances between sequences, which is calculated by `from_fasta_to_distance_matrix.py`.

The `updated_cactus_code_for_hiv_data.py` has been specialized for visualization of the result, users can apply it to `hiv_seq_hamming_distance.csv`

### For 20 random points and color circle data

The generate method for 20 random points is shown in `generate_random_integer_points.py`, and the one of the possible results has been shown in `manhattan_distances_20x20.csv`.

The `updated_cactus_code_for_20x20_points_and_colorcircle.py` has been specialized for visualization of the result, users can apply it to `manhattan_distances_20x20.csv` and `color_matrix_output.csv`
