#!/bin/bash

cluster_id="$1"
ipengine --cluster-id=$cluster_id &
sleep 10
python test_ipyparallel.py
