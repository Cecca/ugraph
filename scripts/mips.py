#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
from pprint import pprint
import itertools
import networkx as nx


def remap_protein(protein):
    return protein.upper()


def load_proteins(path):
    proteins = set()
    with open(path, "r") as fp:
        for line in fp.readlines():
            u, v, prob = tuple(line.split())
            proteins.add(remap_protein(u))
            proteins.add(remap_protein(v))
    return proteins


def load_mips(path):
    complexes = dict()
    with open(path, "r") as fp:
        for line in fp.readlines():
            tokens = line.split("|")
            protein = tokens[0]
            cmplx = tokens[1]
            if cmplx not in complexes:
                complexes[cmplx] = set()
            complexes[cmplx].add(remap_protein(protein))
    return complexes


def print_stats(mips_data):
    sizes = np.array([len(c) for c in mips_data.values()])
    print("\tMin size   ", np.min(sizes))
    print("\tMedian size", np.median(sizes))
    print("\tMean size  ", np.mean(sizes))
    print("\tMax size   ", np.max(sizes))


def filter_complexes(mips_data, proteins):
    new_data = dict()
    for cmplx in mips_data:
        new_complex = mips_data[cmplx].intersection(proteins)
        if len(new_complex) > 0:
            new_data[cmplx] = new_complex
    return new_data


def build_pairs(mips_data):
    pairs = set()
    num_pairs = 0
    counted_pairs = 0
    for cmplx in mips_data:
        cmplx_proteins = mips_data[cmplx]
        num_pairs += (len(cmplx_proteins) * (len(cmplx_proteins) - 1)) / 2
        for pair in itertools.combinations(cmplx_proteins, 2):
            pairs.add(tuple(sorted(pair)))
            counted_pairs += 1
    print(num_pairs, counted_pairs, len(pairs))
    return pairs


def main():
    argp = argparse.ArgumentParser("mips.py")
    argp.add_argument(
        "--ppi",
        help="Protein-Protein interaction network",
        metavar="FILE",
        required=True)
    argp.add_argument(
        "--output-ground", help="Output ground file", metavar="FILE", required=True)
    argp.add_argument(
        "--output-graph", help="Output graph file", metavar="FILE", required=True)
    argp.add_argument(
        "--mips",
        help="MIPS complexes database",
        metavar="FILE",
        required=True)
    args = argp.parse_args()

    proteins = load_proteins(args.ppi)
    print("Loaded", len(proteins), "proteins")
    mips_data = load_mips(args.mips)
    print("Complexes before filtering are", len(mips_data))
    print_stats(mips_data)
    mips_filtered = filter_complexes(mips_data, proteins)
    print("Complexes after filtering are", len(mips_filtered))
    print_stats(mips_filtered)
    pairs = build_pairs(mips_filtered)
    print("There are", len(pairs), "interacting pairs")
    ground_proteins = set()
    for u, v in pairs:
        ground_proteins.add(u)
        ground_proteins.add(v)
    valid_proteins = proteins.intersection(ground_proteins)
    print("Intersection", len(valid_proteins))
    with open(args.output_ground, "w") as fp:
        for p in pairs:
            fp.write("{} {}\n".format(p[0], p[1]))
    G = nx.read_edgelist(args.ppi, data=[('probability', float)])
    G = G.subgraph(valid_proteins)
    with open(args.output_graph, "w") as fp:
        for u, v, p in G.edges_iter(data='probability'):
            fp.write("{} {} {}\n".format(u, v, p))
    


if __name__ == '__main__':
    main()
