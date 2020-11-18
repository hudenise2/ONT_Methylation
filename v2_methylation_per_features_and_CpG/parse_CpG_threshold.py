#!/usr/bin/python3

import argparse, re
import sys, logging

logging.basicConfig(filename='parse_CPG10_olog', level=logging.INFO, format='%(asctime)s %(message)s')
def get_average_per_window(prob_pos_dic):
    new_dic={}
    for k in sorted(prob_pos_dic.keys()):
        new_k=transform_coordinate(str(k))
        if new_k in new_dic:
            new_dic[new_k].append(int(prob_pos_dic[k]))
        else:
            new_dic[new_k]=[int(prob_pos_dic[k])]
    return_dic={k:(int(sum(v)/len(v))) for k,v in new_dic.items()}
    return return_dic

def open_fastq_file(fastq_file):
    open_file=open(fastq_file, 'r')
    return_dic={}
    count=0
    for line in open_file:
        if line.startswith("@") and "-" in line:
            seq_name=line.rstrip()[1:]
            count=0
        elif count==0:
            return_dic[seq_name]={}
            #for dn in ('CG', 'CC', 'CA', 'CT', 'TG', 'GG', 'AG'):
            pos_list=[p.start() for p in re.finditer("CG", line)]
            return_dic[seq_name]["CG"]=pos_list
            count+=1
    return return_dic

def transform_coordinate(coordinate):
    new_coordinate=""
    sign=""
    if coordinate!= 'begin':
        if abs(int(coordinate)) != int(coordinate):
            sign=coordinate[0]
            coordinate=coordinate[1:]
        if len(coordinate) >2:
            if int(coordinate[-2:]) >= 50:
                new_coordinate=sign+str(int(coordinate[:-2])+1)+'01'
            else:
                new_coordinate=sign+coordinate[:-2]+'01'
        else:
            new_coordinate=sign+'01'
    return new_coordinate

def parse_mod(start_dic):
    mod_pos_dic={}
    flag=0
    for line in open("Table_mod", "r"):
       line=line.rstrip()
       dataline=line.split()
       if len(dataline)==1:
           seq=dataline[0]
           countline=0
           flag=0
           if seq in start_dic:
               mod_pos_dic[seq]={}
               flag=1
               logging.info(seq)
       else:
           countline+=1
           if flag==1 and int(dataline[3]) > 192:
               mod_pos_dic[seq][countline+start_dic[seq]]=int(dataline[3])
    return mod_pos_dic

def filter_mod_dic(mod_pos_dic, start_dic):
   dint_dic={}
   filtered_mod_pos_dic={}
   for seq in mod_pos_dic:
       #for each reads in reads: get position of dinucleotide
       fastq_file=fastq_dic[seq]
       logging.info("filtering seq "+seq+" found in fastq file "+fastq_file)
       if fastq_file not in dint_dic:
          dint_dic[fastq_file]=open_fastq_file(fastq_file)
       #for each reads: get the CG position
       logging.info("   get CG position")
       list_CG=list(x+start_dic[seq] for x in dint_dic[fastq_file][seq]['CG'])
       logging.info("   filtering")
       filtered_mod_pos_dic[seq]={x:y for x,y in mod_pos_dic[seq].items() if x in list_CG}
   return filtered_mod_pos_dic

def parse_file(filename):
    cpg_index=-1
    #dictionary of read associated with a given CpG
    read_in_CPG={}
    #dictionary of coordinates of each CpG
    CPG_dic={}
    #dictionary of coordinates of the mapped region for each read
    read_CpGdic={}
    #dictionary of coordinates of the unmapped region for each read
    read_offdic={}
    #dictionary of start position of each reads
    read_pos={}
    for line in open(filename, "r"):
        line=line.rstrip()
        #get the coordinates of the CpG islands
        if line.startswith('chr'):
            cpg_index+=1
            coordinates=line.split(":")[1].split()[0].split("-")
            if coordinates[0] != 'begin':
                new_coordinates=[int(transform_coordinate(coordinate)) for coordinate in coordinates]
                CPG_dic["CpG"+str(cpg_index)]=new_coordinates
            read_CpGdic["CpG"+str(cpg_index)]={}
            read_offdic["CpG"+str(cpg_index)]={}
        else:
            dataline=line.split()
            read=dataline[0]
            start_match=int(dataline[3])
            if "CpG"+str(cpg_index) not in read_in_CPG:
                read_in_CPG["CpG"+str(cpg_index)]=[]
            read_in_CPG["CpG"+str(cpg_index)].append(read)
            read_coordinates=[str(start_match), str(start_match+int(dataline[4]))]
            new_read_coordinates=[int(transform_coordinate(coordinate)) for coordinate in read_coordinates]
            start_match=new_read_coordinates[0]
            if len(dataline)==5:
                #there are not start or end soft-clipping
                #start of mapping is the start as indicated in bam extraction
                #end of mapping is the closest to the start: end of CPG or end of mapped region
                end_match=int(new_read_coordinates[1])
                if int(new_read_coordinates[1]) > CPG_dic["CpG"+str(cpg_index)][1]:
                    end_match= CPG_dic["CpG"+str(cpg_index)][1]
                unmapped_end_read=[end_match, int(transform_coordinate(str(new_read_coordinates[1])))]
                mapped_read=[start_match, end_match]
                if read not in read_CpGdic["CpG"+str(cpg_index)]:
                    read_CpGdic["CpG"+str(cpg_index)][read]=[mapped_read]
                else:
                    if mapped_read not in read_CpGdic["CpG"+str(cpg_index)][read]:
                        read_CpGdic["CpG"+str(cpg_index)][read].append(mapped_read)
                if (unmapped_end_read[1]-unmapped_end_read[0]) !=0:
                    if read not in read_offdic["CpG"+str(cpg_index)]:
                        read_offdic["CpG"+str(cpg_index)][read]=[]
                    if unmapped_end_read not in read_offdic["CpG"+str(cpg_index)][read]:
                            read_offdic["CpG"+str(cpg_index)][read].append(unmapped_end_read)
                read_pos[read]=start_match
            else:
                #there are start and/or end soft-clipping
                start_sclip=0
                end_sclip=0
                read_length=int(dataline[4])
                if "M" in dataline[5]:
                    dataline[5]=dataline[5][dataline[5].index("M")+1:]
                #case where there is start soft-clipping
                if dataline[5] != "-":
                    logging.info(dataline)
                    start_read=new_read_coordinates[0] - int(transform_coordinate(dataline[5][:-1]))+1
                else:
                    start_read=new_read_coordinates[0]
                #case where there is end soft-clipping
                if len(dataline)==7:
                    end_clipping_string=dataline[6]
                    #ensure we only keep the soft-clipping part
                    if "M" in end_clipping_string:
                        end_sclip=int(end_clipping_string[end_clipping_string.index("M")+1:-1])
                    else:
                        end_sclip=int(end_clipping_string[:-1])
                if int(new_read_coordinates[1]) > CPG_dic["CpG"+str(cpg_index)][1]:
                    end_match= CPG_dic["CpG"+str(cpg_index)][1]
                else:
                    end_match=int(transform_coordinate(str(new_read_coordinates[1]-end_sclip)))
                unmapped_start_read=[start_read, new_read_coordinates[0]]
                unmapped_end_read=[end_match, int(transform_coordinate(str(new_read_coordinates[1]-end_sclip)))]
                mapped_read=[start_match, end_match]
                if read not in read_CpGdic["CpG"+str(cpg_index)]:
                    read_CpGdic["CpG"+str(cpg_index)][read]=[mapped_read]
                else:
                    if mapped_read not in read_CpGdic["CpG"+str(cpg_index)][read]:
                        read_CpGdic["CpG"+str(cpg_index)][read].append(mapped_read)
                if (unmapped_start_read[1]-unmapped_start_read[0]) !=0:
                    if read not in read_offdic["CpG"+str(cpg_index)]:
                        read_offdic["CpG"+str(cpg_index)][read]=[]
                    if unmapped_start_read not in read_offdic["CpG"+str(cpg_index)][read]:
                            read_offdic["CpG"+str(cpg_index)][read].append(unmapped_start_read)
                if (unmapped_end_read[1]-unmapped_end_read[0]) !=0:
                    if read not in read_offdic["CpG"+str(cpg_index)]:
                        read_offdic["CpG"+str(cpg_index)][read]=[]
                    if unmapped_end_read not in read_offdic["CpG"+str(cpg_index)][read]:
                            read_offdic["CpG"+str(cpg_index)][read].append(unmapped_end_read)
                read_pos[read]=start_read
    #print(CPG_dic)
    #  {'CpG0': [14801, 15601], 'CpG1': [18101, 18701]}
    #print(read_CpGdic)
    #  {'CpG0': {'36e80244-453d-44bf-abb0-3b1371231924': [[14001, 15601]],
    #         '1a3f737a-5743-47a3-9382-c1366a50c025': [[15101, 15601]],
    #         '8cc2f154-d589-4f55-ba82-44be31e4d392': [[15101, 15601]]},
    #  'CpG1': {'36e80244-453d-44bf-abb0-3b1371231924': [[201, 10801]],
    #          '1a3f737a-5743-47a3-9382-c1366a50c025': [[15101, 18701]],
    #          '8cc2f154-d589-4f55-ba82-44be31e4d392': [[15101, 18701]]}}
    #print(read_offfic)
    #  {'CpG0': {'36e80244-453d-44bf-abb0-3b1371231924': [[15601, 25501]],
    #          '1a3f737a-5743-47a3-9382-c1366a50c025': [[12101, 15101], [15601, 29901]],
    #          '8cc2f154-d589-4f55-ba82-44be31e4d392': [[7701, 15101], [15601, 34001]]},
    #  'CpG1': {'1a3f737a-5743-47a3-9382-c1366a50c025': [[12101, 15101], [18701, 30801]],
    #          '8cc2f154-d589-4f55-ba82-44be31e4d392': [[7701, 15101], [18701, 34001]]}}
    #print(read_pos)
    #  {'36e80244-453d-44bf-abb0-3b1371231924': 201,
    #  '1a3f737a-5743-47a3-9382-c1366a50c025': 12101,
    #   '8cc2f154-d589-4f55-ba82-44be31e4d392': 7701}
    return CPG_dic, read_CpGdic, read_offdic, read_pos

filename=sys.argv[1]
# read position of CPG_islands & read position of read_positions (clipping or not) in relation to chr
CPG, reads, clipped_reads, start_dic= parse_file(filename)
logging.info("read chr_CPG_reads file")
f=filename.split("/")[-1].split("_")[0]
#for each reads in reads: get the fastq file it belongs to
fastq_dic={}
for line in open("FASTQ_dic", "r"):
    line=line.rstrip()
    dataline=line.split(",")
    fastq_dic[dataline[1]]=dataline[0]
logging.info("read FASTQ_dic file")

#parse the modified table file
filtered_mod_pos_dic={}
mod_pos_dic=parse_mod(start_dic)
logging.info("read modified_table file")
#filter the mod_pos_dic for CG dinucleotide
filtered_mod_pos_dic=filter_mod_dic(mod_pos_dic, start_dic)
logging.info("parse modified_table file")

#for each CPG island read, get the CG position
CG_dic={}
new_mod_pos_dic={}
end_dic={}
output_file=open("Threshold"+f+"_output_parse_dic", "w")
for island in reads:
    end_dic[island]={}
    for read in reads[island]:
        end_dic[island]["ON"]=[]
        end_dic[island]["OFF"]=[]
        logging.info("   - got CG position on read "+read)
        if read in mod_pos_dic:
            new_mod_pos_dic[read]=get_average_per_window(filtered_mod_pos_dic[read])
            logging.info("   -got average_per window")
            list_out_CPG=[]
            list_in_CPG=list(range(reads[island][read][0][0], reads[island][read][0][1], 100))
            if read in clipped_reads[island]:
                for region in clipped_reads[island][read]:
                    list_out_CPG+=list(range(region[0], region[1], 100))
                logging.info(str(len(list_in_CPG)))
                logging.info(str(len((list_out_CPG))))

                end_dic[island]["ON"]+=[new_mod_pos_dic[read][str(x)] for x in list_in_CPG if str(x) in new_mod_pos_dic[read]]
                end_dic[island]["OFF"]+=[new_mod_pos_dic[read][str(x)] for x in list_out_CPG if str(x) in new_mod_pos_dic[read]]
                logging.info(str(end_dic[island]))
                for switch in end_dic[island]:
                    output_file.write(island+"\t"+switch+"\t"+str(end_dic[island][switch])+"\n")
logging.info("end of run")
