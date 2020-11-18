#!/usr/bin/python3

#/rds/project/rd109/rds-rd109-durbin-group/projects/cichlid/ONT_Basecall/Hubert_OtoArg_fin_L12_2020-03-16/HDPython3/bin/python

#/usr/bin/python3

import argparse, re, sys, logging 
from datetime import datetime
import statistics as st
import pysam
import random
from Bio import SeqIO


logging.basicConfig(filename='proc_e',  level=logging.DEBUG)

def open_fastq(fastq):
    pos_dic={}
    print(datetime.now())
    with open(fastq) as f:
       for record in SeqIO.parse (f, "fastq"):
          seq_name =record.id
          seq=str(record.seq)
          pos=[x.start() for x in re.finditer("CG", seq)]
          pos_dic[seq_name]=pos
    print(datetime.now())
    return pos_dic

def open_modif(prefix, pos_dic):
   mod_dic={}
   print("\n")
   print(datetime.now())
   with open("modif_tables/Table_mod_"+prefix+"_PASS") as m:
      for line in m:
         line=line.rstrip()
         if "-" in line:
            mseq_name=line
            count=-1
            mod_dic[mseq_name]=[]
         else:
            count+=1
            if count in pos_dic[mseq_name]:
               mod_dic[mseq_name].append(int(line.split()[3]))
   print(datetime.now())
   return mod_dic



fastq=sys.argv[1]
prefix=fastq.split("_")[1]
pos_dic=open_fastq(fastq)
mod_dic=open_modif(prefix, pos_dic)
outfile=open("modif_tables/Preprocessed_Table_mod2_"+prefix+"_PASS", "w")
for entry in mod_dic:
   outfile.write(">"+entry+"\n")
   outfile.write(str(mod_dic[entry])[1:-1]+"\n")
outfile.close()
outfile=open("modif_tables/Preprocessed_Table_pos_"+prefix+"_PASS", "w")
for entry in pos_dic:
   outfile.write(">"+entry+"\n")
   outfile.write(str(pos_dic[entry])[1:-1]+"\n")
outfile.close()
   
