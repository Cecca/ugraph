#!/usr/bin/env python

"""Builds a confusion matrix given a clustering, based on a ground truth.

"""


import gzip
import bz2
import json
import argparse
import sys
from itertools import combinations


def _num_distinct_proteins(clusters):
    proteins = set()
    for cluster in clusters:
        for p in cluster:
            proteins.add(p)
    return len(proteins)

def build_pairs(data):
    if type(data) == list:
        clusters = data
        print("Building pairs from", len(clusters), "clusters, comprising",
              _num_distinct_proteins(clusters), "distinct proteins",
              file=sys.stderr)
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
            "Parameter of wrong type, expected either a list of clusters or a map with a 'clustering' key, was\n", type(data))
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
    """Computes the confusion matrix of `actual_pairs` with respect to `ground_pairs`.

    In the computation, only pairs with both elements contained in the ground set
    are considered, so to have meaningful numbers.
    """
    proteins = set()
    for u, v in ground_pairs:
        proteins.add(u)
        proteins.add(v)

    filtered_actual_pairs = []
    for u, v in actual_pairs:
        if u in proteins and v in proteins:
            filtered_actual_pairs.append((u, v))
    actual_pairs = filtered_actual_pairs
    print("Filtered pairs are", len(actual_pairs), file=sys.stderr)

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
    possible_pairs = len(proteins)*(len(proteins)-1)/2
    print("possible pairs are", possible_pairs, file=sys.stderr)
    negatives = possible_pairs - positives

    tn = negatives - fp

    assert positives == tp + fn
    assert negatives == fp + tn

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
        line_tokens = [l.split() for l in fp.readlines()]
    is_pairs = True
    for tokens in line_tokens:
        if len(tokens) != 2:
            is_pairs = False
            break
    if is_pairs:
        print("Ground truth: `pairs` format", file=sys.stderr)
        return line_tokens
    else:
        print("Ground truth: `one cluster per line` format", file=sys.stderr)
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
    if not isinstance(raw, str):
        raw = raw.decode('utf-8')
    if raw.startswith('{'):
        print("Actual: `json` format", file=sys.stderr)
        data = json.loads(raw)
        return build_pairs(data)
    else:
        print("Actual: `plain text` format", file=sys.stderr)
        line_tokens = [set(l.split()) for l in raw.split("\n")]
        is_pairs = True
        for tokens in line_tokens:
            if len(tokens) != 2:
                is_pairs = False
                break
        if is_pairs:
            print("Actual: `one pair per line` format", file=sys.stderr)
            return line_tokens
        else:
            print("Actual: `one cluster per line` format", file=sys.stderr)
            return build_pairs(line_tokens)


def confusion_matrix_paths(actual_path, ground_path):
    print("Loading ground truth", file=sys.stderr)
    ground_pairs = _load_ground(ground_path)
    print("Loaded ground truth with", len(ground_pairs), "pairs", file=sys.stderr)
    
    print("Loading actual pairs", file=sys.stderr)
    actual_pairs = _load_clustering(actual_path)
    print("Loaded {} pairs".format(len(actual_pairs)), file=sys.stderr)
    result = confusion_matrix(actual_pairs, ground_pairs)
    return result


if __name__ == '__main__':
    argp = argparse.ArgumentParser()
    argp.add_argument('--actual', required=False, default=None)
    argp.add_argument('--ground', required=True)
    args = argp.parse_args()
    result = confusion_matrix_paths(args.actual, args.ground)
    json.dump(result, sys.stdout)
    print("")
    
