#!/bin/bash

docker build -t ugraph-repro .
test -d reproducibility-results || mkdir reproducibility-results
docker run \
  --mount type=bind,source="$(pwd)"/reproducibility-results,target=/results \
  ugraph-repro

