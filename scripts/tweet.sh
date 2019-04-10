#!/bin/bash
TWEET=$1
TWEETLEN=${#TWEET}
if [ $TWEETLEN -eq 0 ] || [ $TWEETLEN -gt 140 ]; then
    if [ $TWEETLEN -gt 140 ]; then
        let EXTRA=$TWEETLEN-140
        echo "Usage: tweet \"message\" (140 chars or less, you're $EXTRA over)"
    else
        echo "Usage: tweet \"message\" (140 chars or less)"
    fi
    exit 1
else
    curl -u RubenAtSanger:OAXtwitter#189 -d status="$1" http://twitter.com/statuses/update.xml
fi
exit 0
