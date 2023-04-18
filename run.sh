#!/bin/bash

rm /home/padniuk/Desktop/SSbarAnalysis/rootfiles/tmp_root/*.root

Files=../data_analysis/*.root
nRunning=$((1))
nTreads=$((1))

echo "Start of the running"

for file in $Files
do
    echo "Run number: $nRunning"
    ./main.exe $file >/dev/null &
    nRunning=$(($nRunning+$((1))))

    if [ $(($nTreads)) -eq  $(($1)) ]
    then
        wait
        nTreads=$((0))
    fi
    nTreads=$(($nTreads+$((1))))
done

echo "Running is completed"

hadd /home/padniuk/Desktop/SSbarAnalysis/rootfiles/tmp_root/output.root /home/padniuk/Desktop/SSbarAnalysis/rootfiles/tmp_root/*.root

echo "Files have been merged"

