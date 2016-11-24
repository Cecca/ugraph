#!/usr/bin/env python

import sys
import subprocess
import tempfile
import os
import argparse
import shutil
import datetime
import json
import time


def abc_to_matrix(abc_path):
    mcl_path = "input.mcl"
    tab_path = "mapping.tab"
    subprocess.call(["mcxload", "--stream-mirror", "-abc", abc_path, "-o",
                     mcl_path, "-write-tab", tab_path])
    print(os.getcwd(), "listing", os.listdir())
    return mcl_path, tab_path


def run_mcl(mcl_path, inflation):
    result_name = "result"
    subprocess.call(["mcl", mcl_path, "-I", str(inflation), "--write-limit",
                     "-o", result_name])
    return result_name, result_name + "-limit"


def build_clusters(result_path, tab_path):
    outs = subprocess.check_output(
        ["mcxdump", "-icl", result_path, "-tabr", tab_path],
        universal_newlines=True)
    clusters = []
    for line in outs.split("\n"):
        c = set(line.split())
        if len(c) > 0:
            clusters.append(c)
    return clusters


def get_centers(limit_path, tab_path):
    outs = subprocess.check_output(
        ["mcxdump", "-imx", limit_path, "-tab", tab_path, "--dump-pairs"],
        universal_newlines=True)
    centers = set()
    for line in outs.split("\n"):
        tokens = line.split()
        if len(tokens) == 3:
            u, v, score = tuple(tokens)
            if u == v:
                centers.add(u)
    return centers


def build_result(clusters, centers, abc_path, inflation):
    tags = {
        "algorithm": "mcl",
        "inflation": float(inflation),
        "graph": abc_path
    }
    date = datetime.datetime.isoformat(datetime.datetime.now())
    stats = [{"num clusters": len(clusters)}]
    clustering = []
    for cluster in clusters:
        center = None
        for c in centers:
            if c in cluster:
                center = c
                break
        for v in cluster:
            clustering.append({"label": v, "center label": center})
    jobj = {"date": date,
            "tags": tags,
            "tables": {"clustering": clustering,
                       "stats": stats}}
    return jobj


def check_ugraph_scores():
    try:
        subprocess.call(["ugraph-scores", "--help"],
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL)
    except:
        print("Could not find the ugraph-scores executable, is it in your $PATH?")
        sys.exit(1)


def mcl(abc_path, inflation):
    check_ugraph_scores()
    _abc_path_abs = os.path.abspath(abc_path)
    calldir = os.getcwd()
    workdir = tempfile.mkdtemp()
    os.chdir(workdir)
    print("Working in", os.getcwd())
    mcl_path, tab_path = abc_to_matrix(_abc_path_abs)
    start = time.monotonic()
    result, limit = run_mcl(mcl_path, inflation)
    end = time.monotonic()
    clusters = build_clusters(result, tab_path)
    centers = get_centers(limit, tab_path)
    jobj = build_result(clusters, centers, abc_path, inflation)
    jobj["tables"]["performance"] = [{"time": int(1000*(end - start))}]
    os.chdir(calldir)
    fname = jobj["date"] + ".json"
    with open(fname, "w") as fh:
        json.dump(jobj, fh)
    print("Wrote JSON result to", fname)
    subprocess.call(["ugraph-scores",
                     "--graph", abc_path,
                     "--clustering", fname])
    shutil.rmtree(workdir)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("graph", metavar="G", help="The input graph.", nargs=1)
    ap.add_argument(
        "-I",
        metavar="INFLATION",
        help="Inflation",
        default=2.0,
        type=float,
        dest="inflation")

    args = ap.parse_args()
    mcl(args.graph[0], args.inflation)


if __name__ == '__main__':
    main()
