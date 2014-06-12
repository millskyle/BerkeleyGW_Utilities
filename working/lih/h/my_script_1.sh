#!/bin/bash


PWSCF='/home/kmills/espresso-5.0/bin/pw.x  -nk 1 -nd 1 -nb 1 -nt 1 '

cd 1-scf
ln ../../*.UPF ./
$PWSCF < in > out



