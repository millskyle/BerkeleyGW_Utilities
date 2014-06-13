#!/bin/bash

currentdir="`pwd`"


PREFIX=`tr -d '\n' < prefix`
echo "Prefix set to $PREFIX"

PWSCF='/home/kmills/espresso-5.0/bin/pw.x  -nk 1 -nd 1 -nb 1 -nt 1 '
PW2BGW='/home/kmills/espresso-5.0/bin/pw2bgw.x '



DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"
echo "Root directory is $DIR"

cd 1-scf
ln -s ../../*.UPF ./
echo "Running first PWSCF calculation"
$PWSCF < in > out
cd ..
ln -s "../1-scf/${PREFIX}.save/" "2-wfn/"
cd 2-wfn/
echo "Running second PWSCF calculation"
$PWSCF < in > out
echo "Converting PWSCF binary wavefunction to BerkeleyGW WFN"
$PW2BGW < pp_in > pp_out


mv wfn.complex WFN_$PREFIX



cd $currentdir
