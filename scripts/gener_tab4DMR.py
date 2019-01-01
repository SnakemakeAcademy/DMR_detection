import sys
import pysam
from pysam import VariantFile
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help = "Input list", required=True)
parser.add_argument("-p", "--prefix", help = "Output file", required=True)
parser.add_argument("-f", "--fai", help = "fasta fai", required=True)
parser.add_argument("-w", "--window", type=int, help = "win size", default = 200)
parser.add_argument("-s", "--step", type=int, help = "step size", default = 200)
parser.add_argument("-d", "--depth", help = "Depth cutoff", default=1,type=int)
options = parser.parse_args()

#tbx = VariantFile(options.input)
out_CG = open(options.prefix + ".CG.tab", 'w')
out_CHG = open(options.prefix + ".CHG.tab", 'w')
out_CHH = open(options.prefix + ".CHH.tab", 'w')
tbx = pysam.TabixFile(options.input)
with open(options.fai) as fileIn:
    for line in fileIn:
        line = line.rstrip()
        ele = line.split("\t")
        stt = 1;
        while(stt + options.window < int(ele[1])):
            end = stt + options.window - 1
            meth_site = defaultdict()
            meth_site["CG"] = [0, 0]
            meth_site["CHG"] = [0, 0]
            meth_site["CHH"] = [0, 0]
            for row in tbx.fetch(ele[0], stt, end):
                tem_ele = str(row).split("\t")
                #Chr9    10      +       0       0       CHH     CTT
                #print(row)
                if "N" in  tem_ele[5]:
                    continue
                #print(row)
                meth_site[tem_ele[5]][0] += int(tem_ele[3])
                meth_site[tem_ele[5]][1] += int(tem_ele[4])
            #print("XX\n")
            for context in ["CG", "CHG", "CHH"]:
                #print(context)
                if meth_site[context][0] + meth_site[context][1] < options.depth:
                    continue
               
                meth_lev = meth_site[context][0]*1.0/(meth_site[context][0] + meth_site[context][1])
                out_str = "{}\t{}\t{}\t{}\t{}\t{}\n".format(ele[0], stt, end, meth_site[context][0], meth_site[context][1], meth_lev)
                #print(out_str)
                if context == "CG":
                    out_CG.write(out_str)
                if context == "CHG":
                    out_CHG.write(out_str)
                if context == "CHH":
                    out_CHH.write(out_str)
            #print("{}\t{}\t{}".format(ele[0], stt, end))
            stt = stt +options.step
