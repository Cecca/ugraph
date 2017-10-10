#!/usr/bin/env python

import networkx as nx
import sys
import random


def load(path):
    G = nx.read_edgelist(path, data=[('probability', float)])
    return G


def pKwikCluster(graph):
    clusters = []
    # Filter edges with probability less than 0.5
    valid_edges = [(u, v) for u, v, d in graph.edges_iter(data='probability')
                   if d >= 0.5]
    G = nx.Graph(valid_edges)
    print('Filtered graph has',G.number_of_nodes(),
          'nodes and', G.number_of_edges(), 'edges',
          file=sys.stderr)
    while G.number_of_nodes() > 0:
        root = random.choice(G.nodes())
        star = list(G[root].keys())
        cluster = [root] + star
        clusters.append(cluster)
        G.remove_nodes_from(cluster)
    return clusters


if __name__ == '__main__':
    G = load(sys.argv[1])
    print('loaded graph with',G.number_of_nodes(),
          'nodes and', G.number_of_edges(), 'edges',
          file=sys.stderr)
    clusters = pKwikCluster(G)
    num_singletons = len([c for c in clusters if len(c) == 1])
    print("Found", len(clusters), "clusters, of which",
          num_singletons, "singletons",
          file=sys.stderr)
    clusters = sorted(clusters, key=lambda c: len(c), reverse=True)
    for cluster in clusters:
        print(" ".join(cluster))




