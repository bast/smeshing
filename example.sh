#!/usr/bin/env bash

./smesh --boundary=/home/bast/tmp/smeshing/data/happy-bear/boundary.txt \
        --islands=/home/bast/tmp/smeshing/data/happy-bear/islands.txt \
        --config=/home/bast/tmp/smeshing/data/happy-bear/config.yml \
        --resolution_function=/home/bast/tmp/smeshing/data/happy-bear/resolution_function.py \
        --output=data.txt
