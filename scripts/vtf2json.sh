#!/bin/bash

grep -v "^#" $1 | sed -e "s/^ *//" | tr -d "\n\t"

exit 0
