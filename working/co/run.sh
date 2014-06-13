#!/bin/bash
currentdir="`pwd`"






PWSCF='mpirun -n 2 /home/kmills/espresso-5.0/bin/pw.x  -nk 1 -nd 1 -nb 1 -nt 1 '
PW2BGW='/home/kmills/espresso-5.0/bin/pw2bgw.x '






DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"
echo "Root directory is $DIR"


if [[ -f prefix ]]; then


PREFIX=`tr -d '\n' < prefix`
echo "Prefix set to $PREFIX"



cd $DIR 

for inputfile in `\ls */*in*`; do
   sed "s|!PREFIX!|$PREFIX|g" < $inputfile > ${inputfile}_new
done

cd 1-scf
ln -s ../../*.UPF ./
echo "Running first PWSCF calculation"
$PWSCF -in in_new &> out
cd $DIR
ln -s "../1-scf/${PREFIX}.save/" "2-wfn/"
cd 2-wfn/
ln -s ../../*.UPF ./
echo "Running second PWSCF calculation"
$PWSCF -in in_new &> out
echo "Converting PWSCF binary wavefunction to BerkeleyGW WFN"
$PW2BGW < pp_in_new > pp_out

cd $DIR
cd 2-wfn

mv wfn.complex ../WFN_$PREFIX



cd $currentdir


else 

echo -e "\n\n ERROR: You must create a prefix file.  The prefix file defines the system, and is used for naming.  It doesn't have to coincide with the molecules in the system.  See the README for more information.  You can create this file by doing: 
   echo 'systemname' > prefix
"

fi



