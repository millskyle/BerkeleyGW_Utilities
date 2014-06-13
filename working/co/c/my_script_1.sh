#!/bin/bash


PREFIX='co'


PWSCF='/home/kmills/espresso-5.0/bin/pw.x  -nk 1 -nd 1 -nb 1 -nt 1 '
PW2BGW='/home/kmills/espresso-5.0/bin/pw2bgw.x '



DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"
pwd
cd 1-scf
ln -s ../../*.UPF ./
$PWSCF < in > out
cd ..
ln -s "../1-scf/${PREFIX}.save/" "2-wfn/"
cd 2-wfn/
$PWSCF < in > out
$PW2BGW < pp_in > pp_out



