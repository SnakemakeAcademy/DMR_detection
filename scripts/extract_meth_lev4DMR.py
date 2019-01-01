import sys
import pysam
from pysam import VariantFile
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--dmr", help = "dmr list", required=True)
parser.add_argument("-t", "--tabix", help = "tabix", required=True)
parser.add_argument("-d", "--depth", help = "Input list", default=1, type=int)
parser.add_argument("--methdiff", type=float, help = "Methylation differnce", default=0.1)
parser.add_argument("-c", "--context", help = "Context [CG CHG CHH]", default=1)
options = parser.parse_args()

tax_list = []
tabix_file = options.tabix.split(",")
for tabix in tabix_file:
    tbx = pysam.TabixFile(tabix)
    tax_list.append(tbx)

with open(options.dmr) as fileIn:
    for line in fileIn:
        line = line.rstrip()
        ele = line.split("\t")
        stt = int(ele[1])
        end = int(ele[2])
        lev = []
        for tem_tabix in tax_list:
            c_num, t_num = 0, 0
            for row in tem_tabix.fetch(ele[0], stt, end):
                tem_ele = str(row).split("\t")
                if tem_ele[5] !=options.context or int(tem_ele[3]) + int(tem_ele[4]) < options.depth:
                    continue
                c_num += int(tem_ele[3])
                t_num += int(tem_ele[4])
            tem_lev = c_num*1.0 /(c_num + t_num)
            lev.append(str(tem_lev))
        if  abs (float(lev[1]) - float(lev[0])) < options.methdiff:
            continue
        out_str = "{}\t{}\t{}\t".format(ele[0], stt, end) + "\t".join(lev)
        print out_str
'''
0 chr
1 stt
2 end
4 context
7 C_num
8 CT_num
'''

