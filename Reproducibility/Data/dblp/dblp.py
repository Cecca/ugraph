#!/usr/bin/env python3

import xmltodict
import gzip
import sys
from itertools import combinations
import time
import math
import tempfile
import networkx as nx
import logging

logging.basicConfig(level=logging.INFO)


def compute_probability(cnt):
    return 1 - math.e**(-0.5 * cnt)


class Counter(object):
    def __init__(self):
        self.count = 0
        self.start = time.time()

    def inc(self, *args):
        self.count += 1
        if self.count % 5000 == 0:
            elapsed = time.time() - self.start
            if elapsed < 60:
                elapsed_str = "{:2.2f} s".format(elapsed)
            elif elapsed < 3600:
                elapsed_str = "{:2.2f} m".format(elapsed / 60)
            else:
                elapsed_str = "{:2.2f} h".format(elapsed / 3600)
            logging.info("++> {:10d} articles | {} | {:6.2f} articles/s".format(
                self.count, elapsed_str, self.count / elapsed))
        return True


class GraphBuilder(object):
    def __init__(self, types):
        self.edges = {}
        self.authors = set()
        self.counter = Counter()
        self.types = types
        self._valid_authors = None

    def valid_article(self, article):
        """Returns true if the article is valid"""
        if 'author' not in article or 'year' not in article:
            return False
        if not isinstance(article['author'], list):
            return False

        if int(article['year']) > 2012:
            return False

        if article['@publtype'] == 'informal publication':
            return False
        
        article_type = article['@key'].split('/')[0]
        if self.types is None:
            return True
        if article_type not in self.types:
            return False
        return True

    def add_edge(self, authA, authB):
        # Add authors
        self.authors.add(authA)
        self.authors.add(authB)
        # Add edge
        if authA < authB:
            key = (authA, authB)
        else:
            key = (authB, authA)
        if key not in self.edges:
            self.edges[key] = 1
        else:
            self.edges[key] += 1

    def get_author_name(self, auth):
        if not isinstance(auth, str):
            return auth["#text"]
        else:
            return auth
            
    def process_article(self, path, article):
        article['@key'] = path[1][1]['key']
        article['@publtype'] = path[1][1]['publtype'] if 'publtype' in path[1][1] else None
        self.counter.inc()
        if self.valid_article(article):
            for authA, authB in combinations(article['author'], 2):
                authA = self.get_author_name(authA)
                authB = self.get_author_name(authB)
                try:
                    self.add_edge(authA, authB)
                except TypeError as e:
                    logging.error("A={}, B={}".format(authA, authB))
                    raise e
        return True

    def parse(self, path):
        with gzip.open(path, 'rb') as fp:
            xmltodict.parse(fp,
                            item_depth=2,
                            xml_attribs=True,
                            item_callback=self.process_article)

    def edges_stream(self):
        for ((a, b), cnt) in self.edges.items():
            yield (a, b, cnt)


def process_name(name):
    return name.replace(" ", "_")


def write_stream(stream, fp):
    for tup in stream:
        fp.write(process_name(tup[0]))
        fp.write(" ")
        fp.write(process_name(tup[1]))
        if len(tup) == 3:
            p = compute_probability(tup[2])
            fp.write(" ")
            fp.write(str(p))
        fp.write("\n")


if __name__ == '__main__':
    logging.info('Parsing input')
    types = {'journals'}
    gb = GraphBuilder(types)
    gb.parse(sys.argv[1])

    logging.info("Building graph")
    G = nx.Graph()
    G.add_weighted_edges_from(gb.edges_stream())
    gb = None # Allow garbage collection
    logging.info("Computing largest connected component")
    LCC = max(nx.connected_component_subgraphs(G), key=len)
    logging.info("Writing the edges to stdout")
    write_stream(LCC.edges_iter(data='weight'), sys.stdout)
