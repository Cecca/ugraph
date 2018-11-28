#!/usr/bin/env python3

import matplotlib
import random
import sys
import experiment
import os
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import rc
import matplotlib.pyplot as plt
from pprint import pprint
import re

REPLACEMENTS = {
    'algorithm': {
        'gmm': r"\textsc{gmm}",
        'mcl': r"\textsc{mcl}",
        'k-center': r"\textsc{mcp}",
        'k-median': r"\textsc{acp}"
    },
    'input': {
        re.compile("collins.*"): 'Collins',
        re.compile("krogan.*"): 'Krogan',
        re.compile("gavin.*"): 'Gavin',
        re.compile("dblp.*"): 'DBLP'
    },
    'measure': {
        'p_min': "$p_{min}$",
        'average probability': "$p_{avg}$",
        'inner-avpr': 'inner-AVPR',
        'outer-avpr': 'outer-AVPR'
    }
}

MCL_FILES = "comparison-results/mcl/*.json*"
EXPERIMENTS_FILES = "comparison-results/*.json.bz2"


def normalize_input(table):
    def _normalize(name):
        f = os.path.basename(name).split(".")[0]
        if f.startswith("collins"):
            return "Collins"
        elif f.startswith("krogan"):
            return "Krogan"
        elif f.startswith("gavin"):
            return "Gavin"
        elif f.startswith("dblp"):
            return "DBLP"
        else:
            return f
    table["input"] = [_normalize(d) for d in table["input"]]

def normalize_algorithm_name(name):
    if name == "gmm":
        return r"$\textsc{gmm}$"
    elif name == "mcl":
        return r"$\textsc{mcl}$"
    elif name == "k-center":
        return r"$\textsc{mcp}$"
    elif name == "k-median":
        return r"$\textsc{acp}$"
    else:
        return name

def measure_label(name):
    if name == "p_min":
        return "$p_{min}$"
    elif name == "average probability":
        return "$p_{avg}$"
    elif name == "time":
        return "time (ms)"
    else:
        return name
    
def normalize_algorithm(table):
    table["algorithm"] = [normalize_algorithm_name(d) for d in table["algorithm"]]


@experiment.cached_table("mcl-scores", MCL_FILES)
def load_mcl():
    table = experiment.load_table(MCL_FILES, "scores")
    print("loaded scores")
    perf = experiment.load_table(MCL_FILES, "performance")[["time"]]
    print("loaded time")
    table = table.join(perf)
    table["input"] = table["graph"]
    table["k"] = table["num clusters"]
    normalize_input(table)
    table = table.replace(REPLACEMENTS)
    table = table[["algorithm", "input", "k", "p_min", "average probability", "time", "inner-avpr", "outer-avpr"]]
    return table

@experiment.cached_table("our-scores", EXPERIMENTS_FILES)
def load_ours():
    table = experiment.load_table(EXPERIMENTS_FILES, "scores")
    perf = experiment.load_table(EXPERIMENTS_FILES, "performance")[["time"]]
    table = table.join(perf)
    normalize_input(table)
    table = table.replace(REPLACEMENTS)
    table = table[["algorithm", "input", "k", "p_min", "average probability", "time", "inner-avpr", "outer-avpr"]]
    return table    


@experiment.cached_table("mcl-proteins-scores", "../MCL/proteins-results/*")
def load_mcl_proteins():
    table = experiment.load_table("../MCL/proteins-results/*", "scores")
    perf = experiment.load_table("../MCL/proteins-results/*", "performance")[["time"]]
    table = table.join(perf)
    table["input"] = table["graph"]
    table["k"] = table["num clusters"]
    normalize_input(table)
    table = table[["algorithm", "input", "k", "p_min", "average probability", "time", "inner-avpr", "outer-avpr"]]
    return table


def format_float_point(val):
    if val == 0.0:
        return "$< 10^{-3}$"
    label = "{:.3f}".format(round(val, 3))[1:]
    label = "${}$".format(label)
    return label


def format_float_exp(val):
    label = "{:.3E}".format(round(val, 3))
    label = label.replace("E", "\cdot 10^{")
    label = "${}{}$".format(label, "}")
    return label


def plot(table, measure, algorithms, sharey=True, exp=False, outname=None):
    print(set(table['input'].values))
    datasets = [
        "Collins",
        "Gavin",
        "Krogan",
        "DBLP"
    ]
    sns.set_palette(sns.color_palette('Paired', len(datasets)))
    algorithms = [normalize_algorithm_name(a) for a in algorithms]
    fig, axs = plt.subplots(nrows=1, ncols=len(datasets),
                            #figsize=(6.97522, 2),
                            figsize=(6.97522, 1.4),
                            sharey=sharey)
    first = True
    if sharey:
        max_height = table.groupby(['k', 'algorithm', 'input']).mean()[measure].max()
    for dataset, ax in zip(datasets, axs):
        plotdata = table[table["input"] == dataset]
        if plotdata.empty:
            print("Empty dataframe for", dataset, ", skipping")
            continue
        sns.barplot(data=plotdata,
                    x="k", y=measure, hue="algorithm",
                    hue_order=algorithms,
                    ax=ax,
                    ci=None)
        
        if not sharey:
            max_height = max([p.get_height() for p in ax.patches])
        max_order = 0
        while max_height / (10.0**max_order) >= 1:
            max_order += 1
        max_order -= 1
        for p in ax.patches:
            xpos = p.get_x() + p.get_width() / 2
            height = p.get_height()
            shift = max_height*0.02
            if not np.isnan(height):# height > 0:
                label = "{:.3f}".format(round(height, 3))[1:]
                label = "${}$".format(label)
                if exp:
                    ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))
                    label = "${:.3f}$".format(height / (10**max_order))
                    xpos += 0.01
                else:
                    label = format_float_point(height)
                if height < 4 * max_height / 5:
                    ypos = height + shift
                    ax.text(xpos, ypos, label, ha="center", va="bottom",
                            color="black", fontsize=7, rotation=90)
                else:
                    ypos = height - shift
                    ax.text(xpos, ypos, label, ha="center", va="top",
                            color="white", fontsize=7, rotation=90)

        ax.set_title(dataset)
        ax.set_xlabel("")
        if first:
            ax.set_ylabel(measure_label(measure))
            first = False
        else:
            ax.set_ylabel("")

    # Reset legends
    handles, labels= None, None
    for ax in axs:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend([], [])

    axs[0].set_xlabel('k', labelpad=1)
    axs[0].xaxis.set_label_coords(0, -0.1)        

    l = fig.legend(handles, labels, ncol=len(labels),
                   loc='lower center',
                   bbox_to_anchor=(0.5, -0.06),
                   frameon=False)
    
    sns.despine(left=True)
    plt.tight_layout(pad=0.0, w_pad=0.6, rect=(0, 0.1, 1, 1))
    plt.savefig("imgs/{}".format(outname.replace(" ", "_")))
    

def to_longform(table):
    table = table.set_index(['algorithm', 'input', 'k'])
    table.columns.name = "measure"
    table = table.stack()
    table.name = 'value'
    table = table.reset_index()
    table = table.replace(REPLACEMENTS)
    return table


def add_labels(ax, max_height=1.0, exp=False):
    for p in ax.patches:
        xpos = p.get_x() + p.get_width() / 2
        height = p.get_height()
        shift = max_height*0.02
        r, g, b, a = p.get_fc()
        if not np.isnan(height):# height > 0:
            label = "{:.3f}".format(round(height, 3))[1:]
            label = "${}$".format(label)
            if exp:
                ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))
                label = "${:.3f}$".format(height / (10**max_order))
                xpos += 0.01
            else:
                label = format_float_point(height)
            if height < 4 * max_height / 5:
                ypos = height + shift
                ax.text(xpos, ypos, label, ha="center", va="bottom",
                        color="black", fontsize='x-small', rotation=90)
            else:
                font_color = "white"
                if r + g + b >= 1.5:
                    font_color = "black"
                ypos = height - shift
                ax.text(xpos, ypos, label, ha="center", va="top",
                        color=font_color, fontsize='x-small', rotation=90)



def plot_multi(table, measures, algorithms, sharey=True, sharex=True, exp=False, outname=None):
    def _plt(**kwargs):
        data = kwargs.pop("data")
        sns.barplot(data=data,
                    x="k", y='value', hue="algorithm",
                    hue_order=algorithms,
                    ci=None)
        dataset = data['input'].values[0]
        
        add_labels(plt.gca(), exp=exp)
        sns.despine(ax=plt.gca(), top=True, bottom=False, left=True, right=True)
        
    datasets = [
        "Collins",
        "Gavin",
        "Krogan",
        "DBLP"
    ]
    sns.set_palette(sns.color_palette('Paired', len(datasets)))
    algorithms = [normalize_algorithm_name(a) for a in algorithms]

    # Bring the table to long form
    table = to_longform(table)
    table = table[table['measure'].isin(measures)]

    width = 6.97522 / len(datasets)
    height = 2.5 / len(measures) # orig submission: 2
    aspect = width / height
    
    g = sns.FacetGrid(data=table, row='measure', col='input', size=height, aspect=aspect,
                      col_order=datasets, row_order=measures, margin_titles=True, sharex=False)
    g.map_dataframe(_plt)
    
    handles = g._legend_data.values()
    labels = g._legend_data.keys()
    
    g.fig.legend(handles, labels, ncol=len(datasets), loc='lower center', bbox_to_anchor=(0,0,1,1))
    g.set_titles(row_template='', col_template='{col_name}')
    for i in range(len(measures)):
        name = g.row_names[i]
        g.axes[i, 0].set_ylabel(name)
        g.axes[i, 0].set_xlabel('k')
        g.axes[i, 0].xaxis.set_label_coords(0, -0.1)
        
        # Reset label names
        for child in g.axes[i, -1].get_children():
            if isinstance(child, matplotlib.text.Annotation):
                child.set_text("")

    plt.tight_layout(rect=(0, 0.1, 1, 1), pad=0.0, w_pad=1, h_pad=1)
    
    if outname is None:
        outname = "{}".format("-".join(measures))
    plt.savefig("imgs/{}".format(outname))


def analysis(table):
    M = REPLACEMENTS['measure']
    A = REPLACEMENTS['algorithm']
    if not os.path.isdir("imgs"):
        os.mkdir("imgs")
    table = table[(table["input"] != 'DBLP') | (table["k"].isin([1818, 5274, 15576]))]
    plot_multi(table, [M["p_min"], M["average probability"]], [A["gmm"], A["mcl"], A["k-center"], A["k-median"]], outname="figure-1.pdf")
    plot_multi(table, [M["inner-avpr"], M["outer-avpr"]], [A["gmm"], A["mcl"], A["k-center"], A["k-median"]], outname="figure-2.pdf")
    plot(table, "time", [A["gmm"], A["mcl"], A["k-center"], A["k-median"]],
         outname="figure-3.pdf",
         sharey=False, exp=True)


def scalability_t(table, algorithms):
    pal = sns.color_palette('Paired')
    sns.set_palette([pal[2], pal[1]])
    A = REPLACEMENTS['algorithm']
    table = table.replace(REPLACEMENTS)
    table = table[table["input"] == 'DBLP']
    table = table[table['k'] < 30000]
    table['time'] = table['time'] / 1000.0
    mcl = table[table['algorithm'] == A['mcl']][['k', 'time']]
    algorithms = [A.get(a, a) for a in algorithms]
    table = table[table['algorithm'].isin(algorithms)]
    print(mcl)
    top = mcl['time'].max()
    print(top)
    
    plt.subplots(figsize=(3.48761, 1.5))
    
    if table.empty:
        print("Empty table in scalability plot")
        return
    ax = sns.pointplot(x='k', y='time', hue='algorithm',
                       data=table,
                       ci=None,
                       #palette=[colors[2], colors[1]],
                       scale=0.6,
                       hue_order=algorithms)

    ax.get_yaxis().get_major_formatter().set_powerlimits((0, 3))
    ax.set_ylabel('time (s)')
    ax.set_xlabel('k', labelpad=1)
    ax.xaxis.set_label_coords(0, -0.05)
    
    ax.plot([0, 1, 2], [top for x in [256, 512, 1024]], marker='X', color='red',
            linestyle=' ', markersize=7)
    
    sns.despine(left=True, bottom=True, offset=256, trim=True)
    plt.tight_layout(pad=0.0, w_pad=1, rect=(0, 0, 1, 0.99))
    plt.savefig('imgs/figure-4.pdf')


    
if __name__ == '__main__':
    A = REPLACEMENTS['algorithm']
    params = {
        'savefig.format': 'svg',
        'font.size': 10,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Helvetica'],
        'text.usetex': True,
    }
    sns.set_context("paper", rc=params)
    sns.set_style("whitegrid")
    pprint(sns.plotting_context())
    rc('text', usetex=True)
    rc('savefig', format='pdf')

    if not os.path.isdir("imgs"):
        os.mkdir("imgs")
    table = pd.concat([load_mcl(), load_ours()])
    table = table[table["k"] > 3]
    #normalize_input(table)
    #normalize_algorithm(table)
    table = table.replace(REPLACEMENTS)

    analysis(table)
    rc('legend', frameon=True)
    scalability_t(table, [A['k-center'], A['mcl']])
