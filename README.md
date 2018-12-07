# Uncertain graphs clustering

This repository contains code and data used in the paper
[_Clustering Uncertain Graphs_](http://www.vldb.org/pvldb/vol11/p472-ceccarello.pdf).

**Contents**

* [Building the code](#building-the-code)
* [Data format](#data-format)
* [Reproducing the results of the paper](#reproducing-the-results-of-the-paper)
* [Result files](#result-files)
* [Citing the work](#citing-the-work)

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

**TL;DR**: run the script `run-docker.sh`

The directory `Reproducibility` contains instructions and data to
reproduce the experimental results described in the paper.

To simplify bulding and running the code on several platforms, a `Dockerfile`
configuration for building a [`docker`](https://www.docker.com/) image is
provided. If you already have `docker`, then running the `run-docker.sh` script
will perform the following actions:

- Build a docker image containing the compiled code and all necessary utilities
  to run experiments and plot figures
- Run the docker image, which will execute the script `Reproducibility/run.sh`
- Copy the resulting files in the newly created directory `reproducibility-results`

In the directory `reproducibility-results/imgs` you will find three PDF files,
named `figure-1.pdf`,  `figure-2.pdf`, and `figure-3.pdf`, each corresponding
to a figure in the paper.

The last console output of the execution will be Table 2 of the paper.

The execution of all the tests will take a long time, so run it on a machine
that can stay dedicated to the task for around one day. If you prefer, you can
change the number of times each experiment is run by changing the first couple
of lines of the script `Reproducibility/run.sh`.

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

## Citing the work

If you find this software useful, please acknowledge its usage by citing 

> Matteo Ceccarello, Carlo Fantozzi, Andrea Pietracaprina, Geppino Pucci, Fabio Vandin. 
> _Clustering Uncertain Graphs_.
> PVLDB, 4(11): 472-44, 2017. DOI: 10.1145/3164135.3164143

```
@article{DBLP:journals/pvldb/CeccarelloFPPV17,
  author    = {Matteo Ceccarello and
               Carlo Fantozzi and
               Andrea Pietracaprina and
               Geppino Pucci and
               Fabio Vandin},
  title     = {Clustering Uncertain Graphs},
  journal   = {{PVLDB}},
  volume    = {11},
  number    = {4},
  pages     = {472--484},
  year      = {2017},
  url       = {http://www.vldb.org/pvldb/vol11/p472-ceccarello.pdf},
  biburl    = {https://dblp.org/rec/bib/journals/pvldb/CeccarelloFPPV17},
}
```
