#!/bin/bash

# starts with a number, the a comma, then another number, and so on:
# 
lines="$1"

regexN='^[[:digit:]]*$'

regexNdN='^([[:digit:]]*)-([[:digit:]]*)$'

regexNCN='^([[:digit:]]*):([[:digit:]]+)$'

regexNcN='^[[:digit:]]*(,?[[:digit:]]*)*$'

regexNC='^([[:digit:]]*):$'

if [[ "$lines" =~ $regexNdN ]]; then

    d1="${BASH_REMATCH[1]}"
    d2="${BASH_REMATCH[2]}"

    #echo "First digit: ${d1} Second digit: ${d2}"

    if [ $d1 -gt $d2 ]; then
        printf "Error: %d bigger than %d\n" "$d1" "$d2"
        exit 1
    fi

    head -n $d2 $2 | tail -n +$d1

    
elif [[ "$lines" =~ $regexNCN ]]; then

    d1="${BASH_REMATCH[1]}"
    d2="${BASH_REMATCH[2]}"

    #printf "numbers: %d  %d\n" "$d1" "$d2"
    
    tail -n +$d1 $2 | head -n $d2

elif [[ "$lines" =~ $regexN ]] || [[ "$lines" =~ $regexNcN ]]; then
    
    lines=${1//,/ }

    for line in `echo $lines`; do  # for line in 9 7 3 8
        head -n $line $2 | tail -1
    done
  
elif [[ "$lines" =~ $regexNC ]]; then

    d1="${BASH_REMATCH[1]}"

    #echo number:$d1
    
    tail -n +$d1 $2
    
fi


exit 0



