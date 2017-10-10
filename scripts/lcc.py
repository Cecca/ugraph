#!/usr/bin/env python

import networkx as nx
import sys


def main():
    inpath = sys.argv[1]
    outpath = sys.argv[2]
    G = nx.read_edgelist(inpath, data=[("p", float)])
    LCC = max(nx.connected_component_subgraphs(G), key=len)
    nx.write_edgelist(LCC, outpath, data=['p'])
    

if __name__ == '__main__':
    main()
