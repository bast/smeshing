#!/usr/bin/env bash

./smesh --boundary=/home/bast/mesh/smeshing/data/fiction/boundary.txt \
        --islands='/home/bast/mesh/smeshing/data/fiction/island*.txt' \
        --config=/home/bast/mesh/smeshing/data/fiction/config.yml \
        --output=/home/bast/mesh/smeshing/data.txt
