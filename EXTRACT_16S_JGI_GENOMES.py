import sys
#*/*.gff would be the sys.argv
from operator import itemgetter
from Bio import SeqIO
import re
with open ("r16s_short.txt","a") as short:
    short.write("IMG Genome ID\tLength\n")
for gff in sys.argv[1:]:
    with open(gff, "r") as annotations:
        r16s_list = []
        for gene in annotations:
            gene = gene.strip()
            if gene.endswith("16S"):
                gene_info = gene.split("\t")
                gene_info_split = gene_info[-1].split(";")
                gene_info_split1 = gene_info_split[0].split("=")
                id=gene_info_split1[1]
                length=int(gene_info[4])-int(gene_info[3])
                orientation=gene_info[6]
                r16s_list.append((id, length,orientation))
        if r16s_list:
            ids=max(r16s_list,key=itemgetter(1))
            print(ids)
            with open("r16s_short.txt","a") as short:
                if ids[1] < 1300:
                    short.write("{}\t{}\n".format(gff[:-15],ids[1]))
            with open("{}_16s.fasta".format(gff[:-4]), "w") as results:
                for record in SeqIO.parse("{}genes.fna".format(gff[:-3]), "fasta"):
                    if record.id == ids[0]:
                        if "+" not in record.description:
                            print(record.id)
                            results.write(">{} {}\n{}\n".format(gff[:-15],record.description,record.seq))
                        else:
                            results.write(">{} {}\n{}\n".format(gff[:-15],record.description,record.seq))
        else:
            with open("16S_not_annotated.txt","a") as noano:
                noano.write("{}\n".format(gff[:-15]))
            with open("r16s_short.txt","a") as short:
                short.write("{}\t0\n".format(gff[:-15]))
