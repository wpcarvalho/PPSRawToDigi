#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: .. <module> <H|V|2H|2V>"
    exit 1
fi

config=$1
typ=$2

RP_combinations_V=("100 104 120 124"
                   "100 104 120"
                   "100 104 124"
                   "104 120 124")

RP_combinations_2V=("100 124"
                    "104 124"
                    "120 124")

RP_combinations_H=("102 103 122 123"
                   "102 103 122"
                   "102 103 123"
                   "103 122 123")

RP_combinations_2H=("102 123"
                    "103 123"
                    "122 123")


if [ "$typ" == "V" ]; then
    RP_combinations=("${RP_combinations_V[@]}")
elif [ "$typ" == "H" ]; then
    RP_combinations=("${RP_combinations_H[@]}")
elif [ "$typ" == "2V" ]; then
    RP_combinations=("${RP_combinations_2V[@]}")
elif [ "$typ" == "2H" ]; then
    RP_combinations=("${RP_combinations_2H[@]}")
else
    echo "Usage: .. <module> <H|V>"
    exit 2
fi

events=1000

for RPs in "${RP_combinations[@]}"
do
    echo "RP combination: "${RPs}
    mkdir -p stats_${RPs// /_}
    for count in 1 2;
    do
        echo $count" particles" && \
        cmsRun set_particle_count_cfg.py $config $count $events $RPs &> /dev/null && \
        mv stats_$count stats_${RPs// /_} && \
        mv hists_${count}.root stats_${RPs// /_} && \
        echo $count" particles finished" &
    done
    wait
done
