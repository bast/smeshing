#!/usr/bin/env bash

./smesh --boundary=/home/bast/mesh/smeshing/data/lofoten/boundary.txt \
        --islands='/home/bast/mesh/smeshing/data/lofoten/islands/*.txt' \
        --config=/home/bast/mesh/smeshing/data/lofoten/config.yml \
        --output=/home/bast/mesh/smeshing/data.txt
