#!/bin/bash

rm -r 1-scf/co.save 1-scf/*UPF
rm -r 2-wfn/co.save 2-wfnf/*UPF
rm 1-scf/out 2-wfn/pp_out 2-wfn/out
rm -r 1*/*.save 2*/*.save
rm 2-wfn/*.wfc 2-wfn/*.UPF 2-wfn/*.complex 2-wfn/*.dat