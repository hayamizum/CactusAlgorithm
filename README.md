# Project Introduction

This project aims at computing a weighted graph that (exactly or approximately) realizes a given distance matrix. For a distance matrix $D$ on a set $X$ of $n$ objects, a realization of $D$ is defined to be a weighted graph $G = (V,E; w)$ such that $X$ is a subset of $V$ and for any two elements $x_i$, $x_j$ in $X$, the shortest path distance $d_G$ between $x_i$ and $x_j$ equals to the input pairwise distance $D(x_i,x_j)$. In the case when $D$ can be realized by a cactus graph or a fully labeled graph, the code outputs a (unique) optimal realization of $D$ (i.e., an exact realization that minimizes the sum of the edge-weights).

Here is a brief example showing what this project aims to do:

![image](https://github.com/hayamizum/CactusAlgorithm/assets/31284538/3e6c8d70-be72-4753-9fcc-e9667ef840ae)

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

The code files are in `Code` directory and data mentioned are in `Data` directory.

When running the code, you can provide a distance matrix in either csv or stdin format.

For csv mode, you need to input a distance matrix as a csv format file that satisfies the following requirements. The absolute path or relative path is needed for searching the file. 

- The csv file must have $n+1$ lines, where $n$ is the number of objects (i.e. the size of the distance matrix).  
- The first line of the csv file only states the number $n$ of objects in the first entry. 
- In the remainder $n$ lines of the csv file, each line states the corresponding row of the distance matrix.

### stdin mode

For stdin mode, it is required to input the distances manually. Shown below is an example using another distance matrix. The first parameter is about the number of objects $n$. Then the input distance matrix including $n×n$ elements is required, and the data must use white space and Enter to distinguish elements. 

```
Input n: 5
Input Distance Matrix: 
0 2 1 1 1
2 0 1 1 1
1 1 0 2 2
1 1 2 0 2
1 1 2 2 0
```

Given the above distance matrix, the code `cactus.py` yields an optimal realization that is a $K_{2,3}$ graph where every edge has weight 1. The code also tells whether or not the generated graph is an exact realization of the given distance matrix (`True` or `False`).

![image](https://github.com/hayamizum/CactusAlgorithm/assets/31284538/0fce0d84-52ac-4f34-a727-8cc3686952a8)

```
True
```

### csv mode

We here demonstrate the code using a sample distance matrix. The sample distance matrix has been created as follows. First, we randomly generated 20 grid points in the plane as below using `generate_random_integer_points.py`, and the points' coordinates are all integers for the reason that the float number will cause precision loss when it comes to add `+` operations. Then, we computed their pairwise distances using L1 metric. The resulting distance matrix is found in `manhattan_distances_20x20.csv`. 

<img width=40% alt="Randomly generated 20 grid points" src="https://github.com/keita1126/CactusAlgorithm/assets/31284538/acefae24-8ed6-47cb-b2e9-9a5dd6a1eaec">

Given the above distance matrix as input, the code `cactus.py` computes the adjacency matrix of its optimal realization and draws it using the Kamada–Kawai layout algorithm. 

```
csv or stdin:csv
File Name:manhattan_distances_20x20.csv
```

<img width=65% alt="output graph for the L1 distance matrix of the random 20 grid points" src="https://github.com/hayamizum/CactusAlgorithm/assets/3113385/bfb22f4d-e59d-4a7e-856c-8a4d0430bc27">

The Kamada–Kawai algorithm generally works well, but as seen below, a slight modification (`updated_cactus_code_for_20x20_points_and_colorcircle.py`) yields a better drawing:

<img width=50% alt="modified output graph for the L1 distance matrix of the random 20 grid points" src="https://github.com/keita1126/CactusAlgorithm/assets/31284538/6a0229ff-2cd6-4aeb-a2c4-7a1297d80ea5">

## Biological data

### p-distances of HIV sequences

HIV sequence data in fasta format are found in `hiv-db_gap_strip.fasta`. Using `from_fasta_to_distance_matrix.py`, you can calculate the p-distances (or Hamming distances) between the sequences, which has been saved in `hiv_seq_hamming_distance.csv`. You can visualize the result by using `updated_cactus_code_for_hiv_data.py`. To make the result more friendly to beginners, the figure has been annotated manually and attached with colors.

#### Attention: If users want to apply to other DNA data, ALL sequences in FASTA file should have the same length

<img width=100% alt="Output graph for the Hamming distances of HIV sequences" src="https://github.com/hayamizum/CactusAlgorithm/assets/3113385/097dac1c-9e5f-443e-acde-4953c232f754">
