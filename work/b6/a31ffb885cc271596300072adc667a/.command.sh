#!/bin/bash -ue
salmon quant --threads 4 --libType=U -i salmon_index -1 gut_1.fq -2 gut_2.fq -o gut
