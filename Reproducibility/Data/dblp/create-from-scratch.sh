#!/bin/bash

# This URL points to the latest DBLP dump
URL=http://dblp.dagstuhl.de/xml/dblp.xml.gz
# In the paper we used the DBLP dump dated 2017-03-03, which can be used
# by uncommenting the following line
# URL=http://dblp.dagstuhl.de/xml/release/dblp-2017-03-03.xml.gz

#curl -O http://dblp.dagstuhl.de/xml/release/dblp-2017-03-03.xml.gz 

RAW_NAME=dblp-raw.xml.gz
if [[ ! -f $RAW_NAME ]]
then
  echo "Downloading DBLP dump to $RAW_NAME (may take some time, depending on your internet connection)"
  curl -o $RAW_NAME http://dblp.dagstuhl.de/xml/release/dblp-2017-03-03.xml.gz 
else
  echo "Raw data already downloaded."
fi

./dblp.py $RAW_NAME > dblp.txt
