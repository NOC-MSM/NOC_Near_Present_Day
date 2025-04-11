#!/bin/bash

set -e

interval=1 # N of seconds between memory queries
record_per_process=0
min_memory=0


usage=$"$(basename "$0") [-h] [-t n] [-p] \n

where: \n
    -h show this help text \n
    -t sampling interval in seconds. Must be an integer \n
    -p record memory per user process. \n
    -m entries with memory below this thresold are ignored. Units in kilobites
"

while getopts m:pt:h flag
do
    case "${flag}" in
        t)  interval=${OPTARG};;
        h)  echo -e $usage 1>&2
            exit 1
            ;;
        p) record_per_process=1 ;;
        m) min_memory=${OPTARG}


    esac
done




if [ $record_per_process == 1 ]
then
    echo "Step Memory Pid Name"
else
    echo "  Step      total    used     free    available"
fi

i=1
while [ 1 ]
do
    if [ $record_per_process == 1 ]
    then
    # saves memory usage per process from ps
        ps -F | tail -n +2 | awk "{ if (\$6 > $min_memory) print $i,\$6,\$2,\$11 }"
    else
    # saves system wide memory usage
        #free --kilo | grep Mem | awk "{if (\$3 > $min_memory) print $i,\$3}"
        free --human | grep Mem | awk "{printf(\"%8.8d %8s %8s %8s %8s\n\",$i,\$2,\$3,\$4,\$7)}"
    fi
    sleep $interval

    i=$((i+interval))

done
