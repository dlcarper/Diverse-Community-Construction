import sys
from Bio import SeqIO
import re
import editdistance
import random
from optparse import OptionParser
from itertools import combinations
parser = OptionParser()
parser.add_option("-i", "--input", dest="input",
                  help="alignment file")
parser.add_option("-D", "--EDITVALUE",type="int",dest="ED_value",
                    help="edit distance you want to examine",default=True)
parser.add_option("-C", "--community",dest="Community_input",
                    help="Starting community")
(options, args) = parser.parse_args()
random.seed(10) #seed for producibility
trim_f= 19 #size of foward primer
sequence_dict={} #empty dictionary will eventually contain sequence ids and sequences
dict_ed={} #empty dictionary that will be a nested dictionaty containing sequence  ids, edvalues and what sequence ids belong in that edit value
with open(options.input) as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"): #loop through every sequence
        #trim primers from alignment
        record.seq=record.seq[trim_f:]
        record.seq=record.seq[:-19]
        #create dictionary where keys are record ids and values are sequences
        sequence_dict[record.id]=str(record.seq)
#compare sequences and get editdistance
#loop through each key and value in the above created dictionary
for key,value in sequence_dict.items():
    #loop through each key and value again
    for key1,value1 in sequence_dict.items():
        #Check if the keys are the same
        if key == key1:
            #if the sequence IDS are the same skip
            continue
        #IF keys are different
        else:
            #Calculate edit distance by comparing each sequence to each other
            edvalue=editdistance.eval(value,value1)
            if (key in dict_ed.keys()):#Check if key is in dictionary
                if (edvalue in dict_ed[key].keys()):#Check if edvalue is in dictionary
                    dict_ed[key][edvalue].append(key1)#Append dictionary with key1
                else:#if edvalue is not in dictionary
                    dict_ed[key][edvalue]=[key1]#Add edvalue and key1 to dicionary
            else:# if key is not in dictionary
                dict_ed[key]={}#Add key to dictionary
                dict_ed[key][edvalue]=[key1]#Add edvalue an key to dictionary
#community =[] # empty list but will eventually contain community members


def withoutcommunityinput():
    community = []
    smallest = 500 #number used that is bigger than number of possibilities
    EDnot_in_dict=[]# empty list but will contain members with no values at edit distance
    for key in dict_ed: #loop through keys in the dictionaty made above
        if options.ED_value not in dict_ed[key]:#Check if the key has an editdistance of input
            EDnot_in_dict.append(key) #if not append the list
        else: # if it does have that edit distqance
            if len(dict_ed[key][options.ED_value])<smallest:#Check if the length of the values for the key at edit distance of input is less than the smallest number
                smallest= len(dict_ed[key][options.ED_value])#set smallest to the length of the smallest list
                smallest_sequence = key#Set smallest seqeunce equal to the key
    if not EDnot_in_dict:# check if this list is empty
        community.append(smallest_sequence) # if it is empty append community with smallest sequence
    else: # if not empty
        community.append(random.choice(EDnot_in_dict))# choose random sequence from list to append community
    print(community)
    print(len(EDnot_in_dict))
    return community

def validateCommunity(starter_community):
    community_validity=[]
    distances=list(range(options.ED_value)) #Create a list of numbers below edit distance value
    distances.append(options.ED_value)# append edit distance value to list
    for pair in combinations(starter_community,2):
        if any([pair[0] in dict_ed[pair[1]][options.ED_value] for options.ED_value in distances if options.ED_value in dict_ed[pair[1]]]):
            print("{},{}".format(pair[0],pair[1]))
            community.append(pair[0])
        else:
            continue
    if community_validity:
        print("community is not valid due to the following members:{}".format(community_validity))
    else:
        community=starter_community
    return community

def loopforCommunity(community):
    not_community=[]#list of members to not put in community
    distances=list(range(options.ED_value)) #Create a list of numbers below edit distance value
    distances.append(options.ED_value)# append edit distance value to list
    while True:# while there are memebers to loop through
        smallest=500
        smallest_sequence=""
        EDnot_in_dict=[]
        for key in dict_ed:# loop through keys agains
            member = False # set to false to keep looping
            if key in community:#check if key is already in the communitiy
                continue #if key in community skip
            else:#if key is not in community check if it is the value for the community members
                for key2 in community:# loop through keys in community
                    if any([key in dict_ed[key2][options.ED_value] for options.ED_value in distances if options.ED_value in dict_ed[key2]]):
                        #check if key you are examining is a member of any member of the communities values at the edit distance 0 to input
                        member=True # set member to true
                        if not key in not_community: # if that key is not in not_community list add it
                            not_community.append(key)
                        break# if it is true break this loop
                if not member: #if key is not a member
                    if options.ED_value in dict_ed[key]:# check if edit value is  present for key
                        if len(dict_ed[key][options.ED_value])<smallest:#Check if the length of the values for the key at edit distance of input is less than the smallest number
                            smallest= len(dict_ed[key][options.ED_value])#set smallest to the length of the smallest list
                            smallest_sequence= key#Set smallest seqeunce equal to the key
                    else: # if edit value isnt present in key
                        EDnot_in_dict.append(key) # add to list
        if smallest_sequence or EDnot_in_dict: # if either of these have a value
            if EDnot_in_dict: # if this list is not empty
                community.append(random.choice(EDnot_in_dict)) # choose random member from list
            else: # if it is empty
                community.append(smallest_sequence) # append community with smallest sequence
        else: # if they dont have a value
            break # break because we are out of members

    return community

if options.Community_input:
    starter_community=[]
    with open(options.Community_input) as starter:
        for line in starter:
            line=line.strip()
            starter_community.append(line)
    community= validateCommunity(starter_community)
    print(community)
    community = loopforCommunity(community)
else:
    community = withoutcommunityinput()
    random.seed(10)
    community = loopforCommunity(community)

print(len(community)) # print length of community
strain_info ={} # empty dicionary for strain information
with open("Strain_INFO.txt") as Info: # open strain file
    for i, line in enumerate(Info): # count lines in file
        line=line.strip()
        if (not i == 0): # if line is not the 1st (so the header)
            fields=line.split("\t") # split on tabs
            strain_info[fields[2]]=fields[0] # set the key of dictionary to genomeID and value to strain
with open ('Community_ED{}.txt'.format(options.ED_value),"w+") as outfile: # open a file to write in
    for key, value in strain_info.items(): # loop through keys and values in strain dictionary
        if key in community: # if that key is found in community
            print("{}\t{}\n".format(key, value))
            outfile.write("{}\t{}\n".format(key, value)) # write the genomeID, a tab, then strain information
