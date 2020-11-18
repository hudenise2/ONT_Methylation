#!/usr/bin/python3
import os, sys
import glob, re
import csv
mod_file=sys.argv[1]  #file with probability -in Table_mod directory
fastq_file=sys.argv[2] #fastq file

def open_fastq_file(fastq_file):
    open_file=open(fastq_file, 'r')
    return_dic={}
    count=0
    for line in open_file:
        if line.startswith("@") and "-" in line:
            seq_name=line.rstrip()[1:]
            count=0
            line_dic={}
        elif count==0:
            return_dic[seq_name]={}
            for dn in ('CA', 'CC', 'CG', 'CT', 'TG', 'GG', 'AG'):
                pos_list=[p.start() for p in re.finditer(dn, line)]
                return_dic[seq_name][dn]=pos_list
            count+=1
    print(list(return_dic.keys())[:10])
    #print(str(return_dic[list(return_dic.keys())[9]]))
    return return_dic

def open_mod_file(fastq_dic, mod_file):
    seq_names=list(fastq_dic.keys())
    open_file=open(mod_file, 'r')
    return_dic={}
    prob_dic={}
    prob_list=[]
    c=0
    for line in open_file:
        line=line.rstrip()
        if re.search("-", line):
            c+=1
            if len(prob_list)>0:
                print(str(c)+" --: "+line)
                prob_dic[name]=prob_list
            name=line
            prob_list=[]
        else:
            if name in seq_names:
                prob_list.append(line.split()[2:4])
    if len(prob_list)!=0: prob_dic[name]=prob_list
    output_file=open("modif_extraction_seq_"+fastq_file, "w")
    for seq_name in prob_dic:
            output_file.write(seq_name+"\n")
            return_dic[seq_name]={}
            for dn in fastq_dic[seq_name]:
                return_dic[seq_name][dn]=[]
                if dn in ['CA', 'CC', 'CG', 'CT']:
                    for entry in fastq_dic[seq_name][dn]:
                        if int(prob_dic[seq_name][int(entry)][0]) < 255:
                            return_dic[seq_name][dn].append(prob_dic[seq_name][int(entry)][1])
                    output_file.write(dn+"\t"+str([prob_dic[seq_name][z][1] for z in fastq_dic[seq_name][dn] if z<=len(prob_dic[seq_name]) and int(prob_dic[seq_name][z][0]) != 255])+"\n")
                else:
                    for entry in fastq_dic[seq_name][dn]:
                        if int(prob_dic[seq_name][int(entry)+1][0]) < 255:
                            return_dic[seq_name][dn].append(prob_dic[seq_name][int(entry)+1][1])
                    output_file.write(dn+"\t"+str([prob_dic[seq_name][z+1][1] for z in fastq_dic[seq_name][dn] if z<len(prob_dic[seq_name])-1 and int(prob_dic[seq_name][z+1][0]) != 255])+"\n")
    return return_dic

def generate_first_output(fastq_dic, mod_dic):
    return_count={x:[0,0] for x in ['CA', 'CC', 'CG', 'CT', 'TG', 'GG', 'AG']}
    for seq_name in mod_dic:
        for dn in mod_dic[seq_name]:
            count_all= len(fastq_dic[seq_name][dn])
            count_mod= len(mod_dic[seq_name][dn])
            return_count[dn]=[return_count[dn][0]+count_all, return_count[dn][1]+count_mod]
    return return_count

def generate_second_output(fastq_dic, mod_dic):
    return_count={x:{} for x in [0,64,128,192]}
    for seq_name in mod_dic:
        for dn in mod_dic[seq_name]:
            for threshold in [0,64,128,192]:
                if dn not in return_count[threshold]: return_count[threshold][dn]=[0,0]
                count_all= len(fastq_dic[seq_name][dn])
                count_mod= len([x for x in mod_dic[seq_name][dn] if int(x) >= threshold])
                return_count[threshold][dn]=[return_count[threshold][dn][0]+count_all, return_count[threshold][dn][1]+count_mod]
    return return_count

os.chdir(".")

#open fastq dic and extract position of dinucleotides
fastq_dic=open_fastq_file(fastq_file)
print("fastq file read")
#open mod table file and extract probability for the position of dinucleotides
mod_dic=open_mod_file(fastq_dic, mod_file)
print("modification file processed")

res=generate_first_output(fastq_dic, mod_dic)

output_file=open("modif_extraction2_"+fastq_file, "w")
output_file.write("dinucleotide\ttotal_count\tmodified_count\tratio_mod/total\tsequences_considered\n")
number_seq=len(mod_dic.keys())
for dn in res:
    output_file.write(dn+"\t"+str(res[dn][0])+"\t"+str(res[dn][1])+"\t"+str(round(res[dn][1]/res[dn][0],2))+"\t"+str(number_seq)+"\n")
output_file.close()
print("first output done")

res2=generate_second_output(fastq_dic, mod_dic)
output_file2=open("modif_extraction_threshold"+fastq_file, "w")
output_file2.write("dinucleotide\tthreshold\ttotal_count\tmodified_count\tratio_mod/total\tsequences_considered\n")
number_seq=len(mod_dic.keys())
for threshold in res2:
    for dn in res2[threshold]:
        output_file2.write(str(threshold)+"\t"+dn+"\t"+str(res2[threshold][dn][0])+"\t"+str(res2[threshold][dn][1])+"\t"+str(round(res2[threshold][dn][1]/res2[threshold][dn][0]*100,2))+"\t"+str(number_seq)+"\n")
output_file2.close()
print("second output done")
