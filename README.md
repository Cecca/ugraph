# Uncertain graphs clustering

This repository contains code for the paper
[_Clustering Uncertain Graphs_](https://arxiv.org/abs/1612.06675).
Specifically, the code in this repository implements
[Algorithm 4](https://arxiv.org/pdf/1612.06675v1.pdf#algocf.4) and
[Algorithm 6](https://arxiv.org/pdf/1612.06675v1.pdf#algocf.6).

## Building the code

To build the code you need

 - CMake >= 3.0
 - Make
 - [peru](https://github.com/buildinspace/peru) to install some
   dependencies
 - OpenMP 
 - [Boost](http://www.boost.org/) >= 1.58
   - graph
   - system
   - filesystem
   - date_time
   - program_options

Then, fromthe root directory of the repository, the following sequence
of commands will build the software

```bash
peru sync                             # get the dependencies
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..   # configure with all optimizations enabled 
make                                  # actually compile the code
```

At this point you should have the following binaries

 - `build/core/ugraph-sc`: implementation of the sequential cluster
   growing strategy of [Algorithm 4](https://arxiv.org/pdf/1612.06675v1.pdf#algocf.4)
 - `build/core/ugraph-cc`: implementation of the concurrent cluster
   growing strategy of [Algorithm 6](https://arxiv.org/pdf/1612.06675v1.pdf#algocf.6)

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

## Usage

### Sequential clustering algorithm

Issuing `ugraph-sc --help` prints a help message to the console. A
typical invocation looks like the following:

```
ugraph-sc --graph path/to/graph --target 256
```

Where `--graph` specifies the path to the graph file, and `--target`
is the number of desired clusters.

### Parallel clustering algorithm

The command `ugraph-cc --help` prints a help message to the
console. Typically, this command is invoked as follows

```
ugraph-cc --graph path/to/graph --batch 256
```

where `--graph` allows to set the input file, and `--batch` is the
number of centers added (on average) in each phase.

### Result files

Both programs produce as a result a json file containing a JSON object
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
