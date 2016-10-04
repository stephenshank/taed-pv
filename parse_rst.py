# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 14:46:49 2016

@author: sshank
"""

import re
from argparse import ArgumentParser

parser = ArgumentParser()
help_string = 'File containing PAML subtree file.'
parser.add_argument('-f', '--file', metavar='FILE', help=help_string, dest='filename')
args = parser.parse_args()
filename = args.filename

p = re.compile('\s+\d+\s+\d+\s+[ATCG-]+.+')
with open(filename, 'r') as input_file:
    for line in input_file:
        if p.match(line):
            print(line, end='')
        if line[:3] == 'Sum':
            break
