#!/bin/bash -ue
salmon quant --threads 4 --libType=U -i salmon_index -1 aml_1.fq -2 aml_2.fq -o aml
