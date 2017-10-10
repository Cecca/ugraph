"""This module is used to handle experimental results"""

import os
from glob import glob
import numpy as np
import json
import dateutil
import pandas as pd
import gzip
import bz2
import matplotlib.pyplot as plt
import seaborn as sns
from functools import wraps


def load_table(globpath, tablename):
    """Load a table from a set of results file"""
    rows = []
    for path in glob(globpath):
        if path.endswith(".gz"):
            fh = gzip.open(path, "rb")
        elif path.endswith(".bz2"):
            fh = bz2.open(path, "rb")
        else:
            fh = open(path, "r")

        try:
            raw = fh.read()
            if isinstance(raw, str):
                json_str = raw
            else:
                json_str = raw.decode("utf-8")
            data = json.loads(json_str)
        except json.JSONDecodeError as e:
            print("Error processing file", path)
            raise e
        except Exception as e:
            print("Error processing file", path)
        if tablename not in data["tables"]:
            # Skip files without the required table
            continue
        try:
            date = dateutil.parser.parse(data['date'])
            tags = data["tags"]
            for data_row in data["tables"][tablename]:
                row = dict()
                row["date"] = date
                row.update(data_row)
                row.update(tags)
                rows.append(row)
        except:
            print("*ERROR* while handling file", path)
        finally:
            fh.close()

    df = pd.DataFrame(rows)
    if len(df) == 0:
        raise Exception("Empty dataframe!")
    df.set_index("date", inplace=True)
    return df


def cached_table(name, input_glob):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            cache_file = ".{}.cache.mpk".format(name)
            if os.path.isfile(cache_file):
                modified_p = False
                mtime = os.path.getmtime(cache_file)
                for f in glob(input_glob):
                    if os.path.getmtime(f) > mtime:
                        modified_p = True
                        break
            else:
                modified_p = True
            if modified_p:
                print("Inputs are newer than the cache file", cache_file)
                table = func(*args, **kwargs)
                if not isinstance(table, (pd.DataFrame, pd.Series)):
                    raise TypeError("Function should return a DataFrame of Series")
                table.to_msgpack(cache_file)
                return table
            else:
                print("Found cached table in file", cache_file)
                table = pd.read_msgpack(cache_file)
                return table
        return wrapper
    return decorator


