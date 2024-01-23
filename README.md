# Project Introduction

This project aims at graph realization of a given two-dimensional distance matrix. A distance matrix contains the distances between all object pairs. For example, if a distance matrix of 10 objects is given, then there should be 10×10 elements. A realization is a graph G(V,E) whose shortest path distances between all pairs of objects are the same as given in the distance matrix. This project focuses on building a minimum realization with the minimum sum of all edges in the graph.

This project is written in Python and can be run in a Python environment.

## Environment Setup

First, clone this repository to your local machine and access the main directory using the command below:

```bash
git clone https://github.com/keita1126/CactusAlgorithm.git
cd CactusAlgorithm
```

To run this project, use:

```bash
python cactus.py
```

and follow the instructions in Input/Output.

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

## Inputs/Outputs

When running the code, there are two modes for input dataset: csv or stdin.

For csv mode, it is required to provide the dataset by csv format file. A path of absolute path or relative path is needed for searching the file.

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
