import sys
from operator import itemgetter
from Bio import SeqIO
import re
import pandas as pd
#sys.argv[1]=Merged_16S.fasta
#primers=515f 806R (work on expanding)
count=0
outfile=open("16S_region.fasta","w")
with open ("Primers_not_present.txt","a") as outfile2:
    outfile2.write("IMG Genome ID\tFW Primer\tRV Primer\tLength\n")
for sequence_file in SeqIO.parse(sys.argv[1],"fasta"):
    forward_matches=re.search(r"GTGCCAGC(A|C)GCCGCGGTAA",str(sequence_file.seq))
    reverse_matches=re.search(r"TTAGA(A|T)ACCC(C|G|T)(A|G|T)GTAGTCC",str(sequence_file.seq))
    forward_matches1=re.search(r"GTGCCAGC(A|C)GCCGCGGTAA",str(sequence_file.seq.reverse_complement()))
    reverse_matches1=re.search(r"TTAGA(A|T)ACCC(C|G|T)(A|G|T)GTAGTCC",str(sequence_file.seq.reverse_complement()))
    if forward_matches and reverse_matches:
        outfile.write(">{} {}\n{}\n".format(sequence_file.id,sequence_file.description,sequence_file.seq[forward_matches.start():reverse_matches.end()]))
        count+=1
    elif forward_matches1 and reverse_matches1:
        sequence_file.seq=sequence_file.seq.reverse_complement()
        outfile.write(">{} {}\n{}\n".format(sequence_file.id,sequence_file.description,sequence_file.seq[forward_matches1.start():reverse_matches1.end()]))
        count+=1
    else:
        if forward_matches:
            with open ("Primers_not_present.txt","a") as outfile2:
                outfile2.write("{}\t{}\tNone\t{}\n".format(sequence_file.id,forward_matches.start(),len(sequence_file.seq)))
        elif reverse_matches:
            with open ("Primers_not_present.txt","a") as outfile2:
                outfile2.write("{}\tNone\t{}\t{}\n".format(sequence_file.id,reverse_matches.end(),len(sequence_file.seq)))
        else:
            with open ("Primers_not_present.txt","a") as outfile2:
                outfile2.write("{}\tNone\tNone\t{}\n".format(sequence_file.id,len(sequence_file.seq)))
print(count)

Info= pd.read_csv("Strain_INFO.txt",sep="\t")
Primers=pd.read_csv("Primers_not_present.txt",sep="\t")
result=pd.merge(Primers,Info, on='IMG Genome ID')
result.to_csv("Needs_resequencing.txt",sep="\t",index=False)
