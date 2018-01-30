#!/usr/bin/env bash

./smesh --boundary=$PWD/data/happy-bear/boundary.txt \
        --islands=$PWD/data/happy-bear/islands.txt \
        --config=$PWD/data/happy-bear/config.yml \
        --output=data.txt
