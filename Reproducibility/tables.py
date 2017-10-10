#!/usr/bin/env python

import numpy as np
import random
import sys
import experiment
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sys.path.append("../scripts")
import confusion
import json
from glob import glob


GROUND, PROTEINS = confusion.load_ground(
    "Data/proteins/TAP_ground.txt")

def confusion_info(group):
    clustering = []
    for idx, row in group[['label', 'center label']].iterrows():
        clustering.append({
            'label': row['label'],
            'center label': row['center label']
        })
    pairs = confusion.build_pairs({'clustering': clustering})
    cm = confusion.confusion_matrix(pairs, GROUND, PROTEINS)
    df = pd.DataFrame([cm])
    return df


@experiment.cached_table("confusion-matrix", "prediction-results/*.json.bz2")
def load():
    table = experiment.load_table("prediction-results/*.json.bz2", "clustering")
    table = table.reset_index().groupby(
        ['date', 'algorithm', 'depth']).apply(confusion_info)
    return table


def load_baseline_krogan():
    with open('prediction-results/krogan-baseline.json') as fp:
        data = json.load(fp)
    data['algorithm'] = 'reference'
    df = pd.DataFrame([data])
    return df

def load_baseline_pkwik():
    rows = []
    with open("prediction-results/pkwik-baseline.json") as fp:
        for line in fp.readlines():
            data = json.loads(line)
            data['algorithm'] = 'pkwik'
            rows.append(data)
    df = pd.DataFrame(rows)
    return df

def load_baseline():
    return pd.concat([load_baseline_krogan(), load_baseline_pkwik()])


def main():
    baseline = load_baseline().groupby("algorithm").mean()[['tpr', 'fpr']]
    print("=== Baseline =========")
    print(baseline)
    #baseline.to_html("tables/baseline.html", border=0, classes=['booktabs'])
    table = load().groupby(level=[1, 2]).mean().rename(index={
        'k-center': 'mcp',
        'k-median-2': 'acp'
    })
    print("=== Our algorithms ===")
    table = table[['tpr', 'fpr']].unstack(0)
    print(table)
    #table.to_html("tables/ours.html", border=0, classes=['booktabs'])



if __name__ == '__main__':
    main()

