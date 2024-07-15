# Optimum Search Schemes for Approximate String Matching Using Bidirectional FM-Index

Welcome to the repository for the paper *"Optimum Search Schemes for Approximate String Matching Using Bidirectional FM-Index"*. This repository contains the implementation and resources related to our research.

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Computer Science Concepts](#computer-science-concepts)
- [Programming Details](#programming-details)
- [Applied Mathematics](#applied-mathematics)
- [Results](#results)
- [Contributors](#contributors)
- [License](#license)

## Introduction

Finding approximate occurrences of a pattern in a text using a full-text index is a central problem in bioinformatics. This repository provides the source code and data for our novel approach to solving this problem using bidirectional FM-index and mixed integer programming (MIP).

## Features

- **Optimal Search Schemes**: Uses MIP to find optimal search schemes for Hamming distance.
- **Bidirectional FM-Index**: Efficiently handles searches from any position within the pattern.
- **Performance**: Significant speed-up over standard backtracking methods.

## More Information

For detailed instructions and usage examples, please refer to the README files inside the following folders:
- `BenchmarkCode-ITV`
- `BenchmarkCode`
- `MIPCode`

## Computer Science Concepts

### Approximate String Matching

The approximate string matching (ASM) problem involves finding substrings in a text that match a pattern approximately, within a specified number of errors (e.g., mismatches, insertions, deletions).

### Bidirectional FM-Index

The bidirectional FM-index allows for flexible searching from any position within the pattern, enabling efficient Approximate String Matching (ASM). It leverages the Burrows-Wheeler Transform (BWT), a key component in many compression algorithms, to support fast substring queries.

The Burrows-Wheeler Transform rearranges the characters of a string into runs of similar characters, which allows for more efficient compression. In the context of string matching, it enables backward searching in the compressed data, facilitating bidirectional searches that start from any position within the pattern and extend in both directions.

The FM-index is built upon the BWT and includes additional auxiliary data structures like the suffix array and rank/select dictionaries. These structures help in performing quick and efficient substring searches:

- **Suffix Array**: An array of integers giving the starting positions of suffixes of a string in lexicographical order. This helps in efficiently finding the position of any substring within the original text.
- **Rank/Select Dictionaries**: These are used to support fast rank and select operations on the BWT. The rank operation counts the number of occurrences of a character up to a certain position, while the select operation finds the position of the nth occurrence of a character.

By combining these elements, the FM-index enables quick navigation through the text in both forward and backward directions. This is particularly useful for approximate string matching, where searches may start from any point within the pattern and require flexibility in handling errors.

### Tries

A search scheme can be visualized by representing each of its searches as a trie that captures all substrings enumerated by the search. Each edge at a level of the trie corresponds to a character of the alphabet at that level of search. A vertical edge represents a match, and a diagonal edge represents a mismatch. The number of edges in the tries provides a measure of the efficiency of the search scheme.

### Dynamic Programming for Edit Distance Calculation

For edit distance calculations, the Smith-Waterman algorithm is employed to handle mismatches, insertions, and deletions efficiently. The dynamic programming approach constructs a matrix where the entry at row i and column j represents the edit distance between the first i characters of the pattern and the first j characters of the text. This matrix can be filled in polynomial time, ensuring that all possible edit distances are considered.

## Programming Details

### Mixed Integer Programming (MIP)

We formulated the search optimization as an MIP problem to find the optimal search schemes. This involves defining variables and constraints to minimize the search steps while covering all possible mismatch patterns.

The MIP formulation is designed to solve the optimal search scheme problem, which is a well-known combinatorial optimization problem. This problem involves finding the search scheme that minimizes the number of steps in ASM-B while ensuring all possible mismatch patterns are covered.

### Key Folders

- `BenchmarkCode-ITV`: Contains the code and scripts for benchmarking the ITV algorithm.
- `BenchmarkCode`: Contains the code and scripts for benchmarking various algorithms.
- `MIPCode`: Contains the implementation of the MIP formulation and solver for finding optimal search schemes.

## Applied Mathematics

### Mathematical Formulation

Our MIP approach is defined by the following formulation:

\[
\text{minimize} \sum_{s=1}^{S} \sum_{l=1}^{R} \sum_{d=0}^{K} n_{s,l,d}
\]

Subject to:

\[
\sum_{i=1}^{P} x_{s,i,j} = 1 \quad \forall s, j
\]

\[
\sum_{j=1}^{P} x_{s,i,j} = 1 \quad \forall s, i
\]

\[
\sum_{i=1}^{P} \sum_{h=1}^{i} x_{s,h,j} - \sum_{i=1}^{P} \sum_{h=1}^{i} x_{s,h,j-1} = t^+_{s,i,j} - t^-_{s,i,j} \quad \forall s, i = 2, \ldots, P-1, j = 1, \ldots, P+1
\]

\[
\sum_{j=1}^{P+1} (t^+_{s,i,j} + t^-_{s,i,j}) = 2 \quad \forall s, i = 2, \ldots, P-1
\]

\[
d - (L_{s,\lceil l/m \rceil} - m\lceil l/m \rceil + l) + 1 \leq (R+1)z_{s,l,d} \quad \forall s, l, d
\]

\[
U_{s,\lceil l/m \rceil} + 1 - d \leq (K+1)z_{s,l,d} \quad \forall s, l, d
\]

\[
\left( \binom{l}{d} (\sigma - 1)^d (z_{s,l,d} + z_{s,l,d} - 2) \right) \leq n_{s,l,d} - n_{s,l-1,d} - (\sigma - 1)n_{s,l-1,d-1} \quad \forall s, l, d
\]

\[
L_{s,i} \leq L_{s,i+1} \quad \forall s, i = 1, \ldots, P-1
\]

\[
U_{s,i} \leq U_{s,i+1} \quad \forall s, i = 1, \ldots, P-1
\]

\[
L_{s,i} + K(\lambda_{q,s} - 1) \leq \sum_{i=1}^{P} \sum_{j=1}^{P} a_{q,j}x_{s,h,j} \leq U_{s,i} + K(1 - \lambda_{q,s}) \quad \forall q, s, i
\]

\[
\sum_{s=1}^{S} \lambda_{q,s} \geq 1 \quad \forall q
\]

\[
n_{s,l,d} \geq 0 \quad \forall q, s, i, j, l, d
\]

\[
L_{s,i}, U_{s,i} \geq 0 \text{ integer} \quad \forall s, i
\]

\[
x_{s,i,j}, \lambda_{q,s}, z_{s,l,d}, t^+_{s,i,j}, t^-_{s,i,j} \in \{0, 1\} \quad \forall q, s, i, j, l, d
\]

## Results

### Performance Gains

Our experiments demonstrate significant improvements over standard backtracking methods, with speed-ups up to 35 times for certain cases. Below are some example results:

![Performance Graph](results/performance.png)

### Example Graphs

Include some visual results from your paper here:

![Optimal Search Scheme](results/optimal_search_scheme.png)

## Contributors

- **Kiavash Kianfar** - Department of Industrial and Systems Engineering, Texas A&M University
- **Christopher Pockrandt** - Department of Computer Science and Mathematics, Freie Universität Berlin
- **Bahman Torkamandi** - Department of Industrial and Systems Engineering, Texas A&M University
- **Haochen Luo** - Department of Industrial and Systems Engineering, Texas A&M University
- **Knut Reinert** - Freie Universität Berlin, Max Planck Institute for Molecular Genetics

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
