#!/bin/bash -ue
salmon quant --threads 4 --libType=U -i salmon_index -1 ecol_1.fq -2 ecol_2.fq -o ecol
