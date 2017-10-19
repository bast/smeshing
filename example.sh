#!/usr/bin/env bash

./smesh --boundary=/home/bast/tmp/smeshing/data/happy-bear/boundary.txt \
        --islands='/home/bast/tmp/smeshing/data/happy-bear/island*.txt' \
        --config=/home/bast/tmp/smeshing/data/happy-bear/config.yml \
        --output=data.txt
