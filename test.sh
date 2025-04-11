#! /bin/bash

for((p=1; p <= 8; p++))
do
    for((k=0; k <= 7; k++))
    do
        for((num=2; num <= 200; num += 10))
        do
            ./a.out -1 1 -1 1 $num 5 $k 1e-14 1000 $p
        done
    done
done