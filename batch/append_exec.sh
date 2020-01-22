#!/bin/bash

eval 'export PATH="/nfs/dust/cms/user/amohamed/anaconda3/bin:$PATH"'

eval 'export KERAS_BACKEND=tensorflow'

source /nfs/dust/cms/user/amohamed/anaconda3/bin/activate hepML

cd $2

/nfs/dust/cms/user/amohamed/anaconda3/envs/hepML/bin/python append_DNN.py --infile $1  --model $3
