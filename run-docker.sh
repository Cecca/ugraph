#!/bin/bash

pushd Reproducibility/Data/dblp
bunzip2 -k dblp-vldb-publication.txt.bz2
popd

docker build -t ugraph-repro .
test -d reproducibility-results || mkdir reproducibility-results
docker run \
  -m 8G \
  --mount type=bind,source="$(pwd)"/reproducibility-results,target=/results \
  ugraph-repro

