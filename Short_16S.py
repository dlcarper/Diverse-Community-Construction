import sys
import pandas as pd
import re

Seq = pd.read_csv("Seq_info.txt",sep="\t")
Info= pd.read_csv("Strain_INFO.txt",sep="\t")
short=pd.read_csv("r16s_short.txt",sep="\t")
result=pd.merge(Info,Seq, on='Strain')
result['Length']=result['16S'].str.len()
results=pd.merge(short,result, on='IMG Genome ID')
results.to_csv("short_16S_withsequence.txt",sep="\t",index=False)

with open ("short_16S_withsequence.txt",'r') as file:
    for i, line in enumerate(file):
        line=line.strip()
        line=line.replace('\"','')
        if (not i == 0):
            fields=line.split("\t")
            fields=list(filter(None,fields))
            if fields[1] <= fields[5]:
                with open("{}/{}_16s.fasta".format(fields[0],fields[0]), 'w') as results:
                    results.write(">{} {}\n{}\n".format(fields[0],fields[2],fields[4]))
