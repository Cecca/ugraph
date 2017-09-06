#!/usr/bin/env python

"""Reconstruct Krogan's et al ground truth.

You can get the file here: http://tap.med.utoronto.ca/exttap/downloads/MIPS_annotations_for_Krogan-etal_Complexes.xls"""

import csv
import sys
from pprint import pprint

def load_data(path):
    with open(path) as fp:
        reader = csv.DictReader(fp)
        return list(reader)


def build_clusters(data):
    clusters_map = dict()
    for row in data:
        cids = row['MIPS complex annotation'].split(",")
        for cid in cids:
            cid = cid.strip()
            if cid != "NOVEL":
                protein = row['ORF']
                if cid not in clusters_map:
                    clusters_map[cid] = {protein}
                else:
                    clusters_map[cid].add(protein)
    return list(sorted(clusters_map.values(), key=len, reverse=True))
    

if __name__ == '__main__':
    path = sys.argv[1]
    outpath = sys.argv[2]
    data = load_data(path)
    clusters = build_clusters(data)
    with open(outpath, "w") as fp:
        for cluster in clusters:
            for p in cluster:
                fp.write(p)    
                fp.write(" ")
            fp.write("\n")

