#!/bin/bash

#set -x

RET=""
SEP=${1:-,}
ENCL=${2:-[]}
OPEN=${ENCL:0:1}
CLOS=${ENCL:1}


while read ENTRY; do
    if [ ! -z "$END_NUM" ]; then
        let NUM=$END_NUM+1
        if [ "$ENTRY" -eq "$NUM" ]; then
            END_NUM=$ENTRY
        else
            if [ "$START_NUM" -ne "$END_NUM" ]; then
                RET+="-${END_NUM}"
            fi
            RET+="${SEP}${ENTRY}"
            END_NUM=$ENTRY
            START_NUM=$END_NUM
        fi
    else
        RET+="${ENTRY}"
        END_NUM=$ENTRY
        START_NUM=$END_NUM
    fi
done < /dev/stdin

if [ "$START_NUM" -ne "$END_NUM" ]; then
    RET+="-${END_NUM}"
fi

RET="${OPEN}${RET}${CLOS}"

echo $RET
