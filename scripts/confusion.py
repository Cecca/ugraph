#!/usr/bin/env python

"""Builds a confusion matrix given a clustering, based on a ground truth.

"""


import gzip
import bz2
import json
import argparse
import sys
from itertools import combinations


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


def confusion_matrix(actual_pairs, ground_pairs):
    tp = 0
    fp = 0
    tn = 0
    fn = 0

    for pair in actual_pairs:
        if pair in ground_pairs:
            tp += 1
        else:
            fp += 1

    for pair in ground_pairs:
        if pair not in actual_pairs:
            fn += 1

    positives = len(ground_pairs)
    proteins = set()
    for u, v in ground_pairs:
        proteins.add(u)
        proteins.add(v)
    possible_pairs = len(proteins)*(len(proteins)-1)/2
    negatives = possible_pairs - positives

    tn = negatives - fp

    return {
        'tpr': tp / positives,
        'fpr': fp / negatives,
        'fnr': fn / positives,
        'tnr': tn / negatives,
        'tp': tp,
        'fp': fp,
        'fn': fn,
        'tn': tn
    }

def _load_ground(path):
    with open(path) as fp:
        line_tokens = [set(l.split()) for l in fp.readlines()]
    is_pairs = True
    for tokens in line_tokens:
        if len(tokens) != 2:
            is_pairs = False
            break
    if is_pairs:
        return line_tokens
    else:
        return build_pairs(line_tokens)


def _clustering_file_handle(path):
    """Get a file handle for gzip files, bz2 files, plain text files, or standard input"""
    if path is None:
        return sys.stdin
    elif path.endswith(".gzip"):
        fp = gzip.open(path)
    elif path.endswith(".bz2"):
        fp = bz2.open(path)
    else:
        fp = open(path)
    return fp


def _load_clustering(path):
    """Load a clustering as a sequence of pairs of nodes belonging to the same cluster.
    
    Works with streamed standard input, (compressed) json files, text files.
    """
    fp = _clustering_file_handle(path)
    raw = fp.read() # Read the whole thing
    if raw[0] == '{':
        data = json.loads(raw)
        return build_pairs(data)
    else:
        line_tokens = [set(l.split()) for l in fp.readlines()]
        is_pairs = True
        for tokens in line_tokens:
            if len(tokens) != 2:
                is_pairs = False
                break
        if is_pairs:
            return line_tokens
        else:
            return build_pairs(line_tokens)


if __name__ == '__main__':
    argp = argparse.ArgumentParser()
    argp.add_argument('--actual', required=False, default=None)
    argp.add_argument('--ground', required=True)
    args = argp.parse_args()
    
    ground_pairs = build_pairs(_load_ground(args.ground))
    
    actual_pairs = _load_clustering(args.actual)
    result = confusion_matrix(actual_pairs, ground_pairs)
    json.dump(result, sys.stdout)
    print("")
