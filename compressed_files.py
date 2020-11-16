# We need to be able to work directly from compressed files (.gz), 
# and to do so without opening them so that they take up less space. 

# We do it like this:

import gzip
with gzip.open('input.gz','r') as infile:
        for line in infile:
            print('got line', line)


# We also should be able to find out whether our input files are compressed
# or not. 

# We do it like this:

import sys
import re

if re.search(r'\.gz', sys.argv[1]):
    with gzip.open('input.gz','r') as infile:
        for line in infile:
            print('got line', line)

else: 
    with open(input.fastq, 'r') as infile:
        for line in infile:
            print('got line', line)


    