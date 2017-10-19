#!/usr/bin/env bash

./smesh --boundary=/home/bast/tmp/smeshing/data/bear/boundary.txt \
        --islands='/home/bast/tmp/smeshing/data/bear/island*.txt' \
        --config=/home/bast/tmp/smeshing/data/bear/config.yml \
        --output=data.txt
