#!/bin/bash

# echo $ALPS_APP_PE

cd $HOME/CA1Sim

cluster_id="test_ipyparallel"
profile='mpi'

set -x

if [ $ALPS_APP_PE == 0 ]; then
    ipcontroller --ip='*' --cluster-id=$cluster_id --profile=$profile --work-dir=$HOME/CA1Sim
elif [ $ALPS_APP_PE == 1 ] ; then
    python test_ipyparallel.py --cluster-id=$cluster_id --profile=$profile --sleep
else
    ipengine --cluster-id=$cluster_id --profile=$profile --work-dir=$HOME/CA1Sim
fi
