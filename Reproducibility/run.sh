#!/bin/bash

## Test python installation
python -c 'import seaborn; import pandas' || {
  echo "Your Python installation is missing some required packages."
  echo "The following are required"
  echo ""
  echo "  seaborn"
  echo "  pandas"
  echo ""
  echo "If you are using the conda distribution, you can setup the exact same"
  echo "python environment used in the paper with the file `conda-environment.yml`"
  echo ""
  echo "Otherwise, just install the missing packages using your distribution's"
  echo "package manager or pip"
  exit 1
}


BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PATH=$PATH:$BASEDIR/../build/core:$BASEDIR/../scripts
NUM_RUNS=1

if [[ ! -f $BASEDIR/Data/dblp/dblp-vldb-publication.txt ]]
then
    echo "Unpacking DBLP dataset"
    bunzip2 $BASEDIR/Data/dblp/dblp-vldb-publication.txt.bz2
fi

##############################################################################
#
# Experiments for Figures 1, 2, 3, and 4
#
##############################################################################

function figures () {
  ## Run comparisons
  COMPARISON_RESULTS_DIR=comparison-results
  test -d $COMPARISON_RESULTS_DIR || mkdir $COMPARISON_RESULTS_DIR
  cd $COMPARISON_RESULTS_DIR
  
  ## Run MCL
  test -d mcl || mkdir mcl
  cd mcl
  
  for INFLATION in 1.2 1.5 2.0
  do
      for GRAPH in $BASEDIR/Data/proteins/collins2007-lcc.txt $BASEDIR/Data/proteins/gavin2006-lcc.txt $BASEDIR/Data/proteins/krogan2006_core-lcc.txt
      do
          mcl.py -I $INFLATION $GRAPH
      done
  done
  
  for INFLATION in 1.15 1.2 1.3
  do
      mcl.py -I $INFLATION $BASEDIR/Data/dblp/dblp-vldb-publication.txt
  done
  
  cd ..
  
  ## Run other algorithms
  
  DATASET=$BASEDIR/Data/proteins/collins2007-lcc.txt
  for TARGET in 24 69 99
  do
      for RUN in $(seq $NUM_RUNS)
      do
          ugraph-gmm --graph $DATASET --target $TARGET --with-avpr
          ugraph-mcpc --graph $DATASET --target $TARGET --rate 0.1 --with-avpr
          ugraph-acpc --graph $DATASET --target $TARGET --rate 0.1 --with-avpr
    done
  done
  
  DATASET=$BASEDIR/Data/proteins/gavin2006-lcc.txt
  for TARGET in 50 172 274
  do
      for RUN in $(seq $NUM_RUNS)
      do
          ugraph-gmm --graph $DATASET --target $TARGET --with-avpr
          ugraph-mcpc --graph $DATASET --target $TARGET --rate 0.1 --with-avpr
          ugraph-acpc --graph $DATASET --target $TARGET --rate 0.1 --with-avpr
      done
  done
  
  DATASET=$BASEDIR/Data/proteins/krogan2006_core-lcc.txt
  for TARGET in 77 289 517
  do
      for RUN in $(seq $NUM_RUNS)
      do
          ugraph-gmm --graph $DATASET --target $TARGET --with-avpr
          ugraph-mcpc --graph $DATASET --target $TARGET --rate 0.1 --with-avpr
          ugraph-acpc --graph $DATASET --target $TARGET --rate 0.1 --with-avpr
      done
  done
  
  DATASET=$BASEDIR/Data/dblp/dblp-vldb-publication.txt
  for TARGET in 77 289 517
  do
      for RUN in $(seq $NUM_RUNS)
      do
          ugraph-gmm --graph $DATASET --target $TARGET --with-avpr
          ugraph-mcpc --graph $DATASET --target $TARGET --rate 0.1 --with-avpr
          ugraph-acpc --graph $DATASET --target $TARGET --rate 0.1 --with-avpr
      done
  done
  
  cd ..

  ./figures.py
}
  
##############################################################################
#
# Experiments for Table 2
#
##############################################################################

function table2 () {
  DATASET=$BASEDIR/Data/proteins/TAP_core.txt
  GROUND=$BASEDIR/Data/proteins/TAP_ground.txt
  
  test -d $BASEDIR/prediction-results/ || mkdir $BASEDIR/prediction-results
  cd $BASEDIR/prediction-results/
  
  # # Evaluate reference clustering
  confusion.py --ground $GROUND --actual $BASEDIR/Data/proteins/TAP_clusters.txt > krogan-baseline.json
  
  # Run pkwik
  test -d pkwik-runs || mkdir pkwik-runs
  for RUN in $(seq 10)
  do
    pkwik.py $DATASET > pkwik-runs/run-$RUN.txt
  done
  for FILE in pkwik-runs/*
  do
    confusion.py --ground $GROUND --actual $FILE >> pkwik-baseline.json
  done
  
  K=547
  RATE=0.1
  for DEPTH in 2 3 4 6 8
  do
    ugraph-mcpc --graph $DATASET --target $K --rate $RATE --depth $DEPTH
    ugraph-acpc --graph $DATASET --target $K --rate $RATE --depth $DEPTH -h 100
  done
  
  cd ..

  ./tables.py
}

## Run everything
figures
table2
