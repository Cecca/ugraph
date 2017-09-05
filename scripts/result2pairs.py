#!/usr/bin/env python

"""Script to convert a json clustering result file to a list of pairs of nodes belonging to the same cluster.

Can be either invoked on the command line, by providing the path to a json file containing the results of a clustering,
or can be used as a library with the method build_pairs.
"""

import gzip
import bz2
import json
import argparse
import sys
from itertools import combinations


def _load_json(path):
    if path.endswith(".gzip"):
        fp = gzip.open(path)
    elif path.endswith(".bz2"):
        fp = bz2.open(path)
    else:
        fp = open(path)

    data = json.load(fp)

    fp.close()
    return data


def build_pairs(data):
    if type(data) == list:
        clusters = data
        pairs = set()
        for c in clusters:
            for u, v in combinations(c, 2):
                pairs.add(tuple(sorted([u, v])))
        return pairs
    if 'clustering' in data:
        table = data['clustering']
    elif 'tables' in data and 'clustering' in data['tables']:
        table = data['tables']['clustering']
    else:
        raise TypeError(
            "Parameter of wrong type, expected either a list of clusters or a map with a 'clustering' key, was\n", data)
    # We have a list of nodes, each with clustering information
    cluster_map = dict()
    for vertex in table:
        center = vertex['center label']
        label = vertex['label']
        if center not in cluster_map:
            cluster_map[center] = [label]
        else:
            cluster_map[center].append(label)
    return build_pairs(list(cluster_map.values()))


if __name__ == '__main__':
    assert len(sys.argv) == 2, "USAGE: result_to_pair.py path.json.bz2"
    path = sys.argv[1]
    data = _load_json(path)
    for u, v in build_pairs(data):
        print(u, v)
