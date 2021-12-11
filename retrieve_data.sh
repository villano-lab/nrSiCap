#!/bin/bash
wget https://osf.io/9xftw/download -O ./data/v3_400k.pkl
wget https://osf.io/xhz5w/download -O ./data/normsi_fast_200k.pkl

# MCMC walks
wget https://osf.io/dw49y/download -O ./data/mcmc_Lind_128walk_50kstep_SNorm_v1.pkl
wget https://osf.io/2nsmh/download -O ./data/mcmc_Lind_128walk_50kstep_SNorm_v2.pkl
wget https://osf.io/25ckn/download -O ./data/mcmc_Sor_128walk_50kstep_SNorm_v1.pkl
wget https://osf.io/a9fnt/download -O ./data/mcmc_Sor_128walk_50kstep_SNorm_v2.pkl
wget https://osf.io/k3xzg/download -O ./data/mcmc_Sor_128walk_50kstep_SNorm_v3.pkl
wget https://osf.io/j7mkg/download -O ./data/mcmc_Sor_128walk_50kstep_SNorm_v4.pkl

# By Series
wget https://osf.io/mgs35/download -O ./data/byseries/07180924_1710/07180924_1710_traces.root
