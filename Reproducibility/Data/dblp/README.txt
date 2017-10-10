= The DBLP dataset =

We use a graph derived from a dump of DBLP, where nodes are authors, and two
authors are connected by an edge if they have co-authored at least one journal
paper together. Each edge has probability 1 - e^(x/2), where x is the number of
collaborations between the two authors.
We then consider the largest connected component of this graph.

In the publication we make use of the graph stored in the file
`dblp-vldb-publication.txt`.
We also provide the script `create-from-scratch.sh` that can be used to
re-create the graph starting from tha raw dump of DBLP.

