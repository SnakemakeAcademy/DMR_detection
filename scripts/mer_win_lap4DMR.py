#!/usr/bin/python

import sys
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help = "Input list", required=True)
parser.add_argument("-o", "--output", help = "Output file", required=True)
parser.add_argument("-d", "--depth", help = "Input list", default=10, type=int)
options = parser.parse_args()

### store windown in dictionary 
dict_reg = defaultdict(dict)
### store table name
list_file = []
file_list = options.input.split(",")
for line in file_list:
    with open(line) as meth_file:
        for methl in meth_file:
            if methl.startswith("#"):
                continue
            methl = methl.rstrip()
            ele = methl.split("\t")
            win_key = '\t'.join(ele[0:3])
            #print '\t'.join(ele[7:])
            if ele[3] + ele[4] < options.depth:
                continue
            c_num = ele[3]
            t_num = ele[4]
            c_sites = "0"  #### Assign 0 
            ele[3] = c_sites
            ele[4] = c_num
            ele[5] = t_num
            if win_key not in dict_reg:
                #print ele[3:] 
                dict_reg[win_key] = ['\t'.join(ele[3:6])]
            else:
                dict_reg[win_key].append('\t'.join(ele[3:6]))

with open(options.output, 'w') as output_file:
    for win_key in dict_reg:
        if len(dict_reg[win_key]) ==  len(file_list): 
            output_file.write(win_key + "\t" + '\t'.join(dict_reg[win_key]) + "\n")
