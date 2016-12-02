#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np
from pprint import pprint
import itertools


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
    for cmplx in mips_data:
        cmplx_proteins = mips_data[cmplx]
        for pair in itertools.combinations(cmplx_proteins, 2):
            pairs.add(tuple(sorted(pair)))
    return pairs


def main():
    argp = argparse.ArgumentParser("mips.py")
    argp.add_argument(
        "--ppi",
        help="Protein-Protein interaction network",
        metavar="FILE",
        required=True)
    argp.add_argument(
        "--output", help="Output file", metavar="FILE", required=True)
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
    with open(args.output, "w") as fp:
        for p in pairs:
            fp.write("{} {}\n".format(p[0], p[1]))
    


if __name__ == '__main__':
    main()
