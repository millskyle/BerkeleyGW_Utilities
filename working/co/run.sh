#!/bin/bash
currentdir="`pwd`"


PREFIX=`tr -d '\n' < prefix`
echo "Prefix set to $PREFIX"

PWSCF='/home/kmills/espresso-5.0/bin/pw.x  -nk 1 -nd 1 -nb 1 -nt 1 '
PW2BGW='/home/kmills/espresso-5.0/bin/pw2bgw.x '






DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"
echo "Root directory is $DIR"



cd $DIR 

for inputfile in `\ls */*in*`; do
   sed "s|!PREFIX!|$PREFIX|g" < $inputfile > ${inputfile}_new
done

cd 1-scf
ln -s ../../*.UPF ./
echo "Running first PWSCF calculation"
$PWSCF < in_new > out
cd $DIR
ln -s "../1-scf/${PREFIX}.save/" "2-wfn/"
cd 2-wfn/
ln -s ../../*.UPF ./
echo "Running second PWSCF calculation"
$PWSCF < in_new > out
echo "Converting PWSCF binary wavefunction to BerkeleyGW WFN"
$PW2BGW < pp_in_new > pp_out

cd $DIR
cd 2-wfn

mv wfn.complex ../WFN_$PREFIX



cd $currentdir
