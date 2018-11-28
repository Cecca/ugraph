# Uncertain graphs clustering

This repository contains code and data used in the paper
_Clustering Uncertain Graphs_.

## Building the code

To build the code you need

 - CMake >= 3.0
 - Make
 - OpenMP
 - [Boost](http://www.boost.org/) >= 1.58 with the following components
   - graph
   - system
   - filesystem
   - date_time
   - program_options

From the root directory of the repository, run the following sequence
of commands to build the software

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..   # configure with all optimizations enabled 
make                                  # actually compile the code
```

At this point you should have the following binaries

 - `build/core/ugraph-mcpc`: implementationof the MCP clustering
   algorithm
 - `build/core/ugraph-acpc`: implementation of the ACP clustering
   algorithm
 - `build/core/ugraph-gmm`: implementation of the classic GMM
   algorithm, adapted to work with probabilities, to be used as a
   baseline.
 - `build/core/ugraph-scores`: executable to compute scores given
   a clustering

## Data format

The software works on undirected graphs, represented in a very simple
textual format. Each line represents an edge, with a source, a
destination, and a probability of existence, as in the following
example

```
0 1 0.64
0 2 0.89
2 4 0.45
1 2 0.13
```

## Reproducing the results of the paper

The directory `Reproducibility` contains instructions and data to
reproduce the experimental results described in the paper.

## Result files

All programs produce as a result a compressed json file containing a JSON object
with the following structure.

 - an object `"tags"` listing the configuration of the run, from the
   git-revision of the code to the input parameters, including the
   seed used for the random number generator
 - the date of the run
 - an object `"tables"` with the following fields:
   - `"scores"`: a list containing a single objects with several
     scores on the quality of the clustering
   - `"clustering"`: a list of objects, where each object represents a
     node of the input graph and has the following fields
     - `"id"`: the internal identifier of the node
     - `"center"`: the internal identifier of the center this node has
       been associated to
     - `"label"`: the label of the node in the input
     - `"center label"`: the label of the center in the input
     - `"probability`": the estimated probability of connection
       between the node and its cluster center
       
These files are suitable to be further processed with software such as [pandas](http://pandas.pydata.org/).

The `scripts/mcl.py` script is a wrapper around the
[MCL](https://micans.org/mcl/) executable that presents the
clusterings computed by MCL in this JSON format, thus allowing uniform
analysis.
