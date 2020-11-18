#!/usr/bin/python3

import argparse, re
import sys

def get_summary_for_all(CPG_dic):
    count_CPG=len(CPG_dic)
    new_list=[0,0,0,0,0,0]
    CPG_islandplus=[]
    CPG_islandminus=[]
    CPG_eq=[]
    for k in CPG_dic:
        for switch in CPG_dic[k]:
            if switch=="ON":
                new_list[0]+=CPG_dic[k][switch][0]
                new_list[1]+=CPG_dic[k][switch][1]
                new_list[2]+=CPG_dic[k][switch][3]
            else:
                new_list[3]+=CPG_dic[k][switch][0]
                new_list[4]+=CPG_dic[k][switch][1]
                new_list[5]+=CPG_dic[k][switch][3]
        if "OFF" in CPG_dic[k]:
            if int(CPG_dic[k]["ON"][2]) == int(CPG_dic[k]["OFF"][2]):
                CPG_eq.append([k, CPG_dic[k]["ON"], CPG_dic[k]["OFF"]])
            elif int(CPG_dic[k]["ON"][2]) > int(CPG_dic[k]["OFF"][2]):
                CPG_islandplus.append([k, CPG_dic[k]["ON"], CPG_dic[k]["OFF"]])
            else:
                CPG_islandminus.append([k, CPG_dic[k]["ON"], CPG_dic[k]["OFF"]])
    new_CPG_islandplus=sort_by_prob(CPG_islandplus)
    new_CPG_islandminus=sort_by_prob(CPG_islandminus)
    new_CPG_islandeq=sort_by_prob(CPG_eq)
    new_list.append(count_CPG)
    new_list.append(len(new_CPG_islandeq))
    new_list.append(len(new_CPG_islandplus))
    new_list.append(len(new_CPG_islandminus))
    output_file=open("CPG/process_CPG_output/"+f+"_INsup_sorted_islands", "w")
    output_file.write("ISLANDS WITH PROBin > PROBout\nisland\tcountCPGin\tprobCmGin\tcountCPGout\tprobCPGout\tdiff\n")
    for entry in new_CPG_islandplus:
        for l in entry:
            output_file.write(l[0]+"\t"+str(l[1][0])+"\t"+str(l[1][2])+"\t"+str(l[2][0])+"\t"+str(l[2][2])+"\t"+str(round(l[1][2] - l[2][2],2))+"\n")
    output_file.close()
    output_file=open("CPG/process_CPG_output/"+f+"_OUTsup_sorted_islands", "w")
    output_file.write("\nISLANDS WITH PROBin < PROBout\n-island\tcountCPGin\tprobCmGin\tcountCPGout\tprobCPGout\tdiff\n")
    for entry in new_CPG_islandminus:
        for l in entry:
            output_file.write(l[0]+"\t"+str(l[1][0])+"\t"+str(l[1][2])+"\t"+str(l[2][0])+"\t"+str(l[2][2])+"\t"+str(round(l[1][2] - l[2][2],2))+"\n")
    output_file.close()
    output_file=open("CPG/process_CPG_output/"+f+"_INeqOUT_sorted_islands", "w")
    output_file.write("\nISLANDS WITH PROBin = PROBout\n-island\tcountCPGin\tprobCmGin\tcountCPGout\tprobCPGout\tdiff\n")
    for entry in new_CPG_islandeq:
        for l in entry:
            output_file.write(l[0]+"\t"+str(l[1][0])+"\t"+str(l[1][2])+"\t"+str(l[2][0])+"\t"+str(l[2][2])+"\t"+str(round(l[1][2] - l[2][2],2))+"\n")
    output_file.close()
    return new_list

def get_median(list):
    if len(list)==1:
        medianval=list[0]
    else:
        sorted_list=sorted(list)
        mid_index=int((len(list)-1)/2)
        if len(list)%2==0:
            medianval=int((sorted_list[mid_index]+sorted_list[mid_index+1])/2)
        else:
            medianval=int(sorted_list[mid_index+1])
    return medianval


def sort_by_prob(CPGlist):
    new_dic={}
    cpg_list=[]
    for entry in CPGlist:
        diff=abs(entry[1][2]-entry[2][2])
        if diff in new_dic:
            if [entry]not in new_dic[diff]:
                new_dic[diff]+=[entry]
        else:
            new_dic[diff]=[entry]
    sort_dic=[new_dic[k] for k in reversed(sorted(new_dic.keys()))]
    return sort_dic

def parse_file(filename):
    cpg_index=-1
    read_xdic={}
    CPG_dic={}
    clip_xdic={}
    sumprob=0
    medianval=0
    for line in open(filename, "r"):
        line=line.rstrip()
        dataline=line.split("\t")
        if dataline[0] not in CPG_dic:
            CPG_dic[dataline[0]]={}
        CPG_dic[dataline[0]][dataline[1]]=[]
        prob=dataline[2][1:-1].split(",")
        CPG_dic[dataline[0]][dataline[1]].append(len(prob))
        sumlist=[int(x.strip()) for x in prob if len(x)>0]
        #if len(sumlist) > 0 : medianval=get_median(sumlist)
        CPG_dic[dataline[0]][dataline[1]].append(sum(sumlist))
        if len(sumlist)==0:
            CPG_dic[dataline[0]][dataline[1]].append(0)
            CPG_dic[dataline[0]][dataline[1]].append(0)
        else:
            CPG_dic[dataline[0]][dataline[1]].append(round(sum(sumlist)/len(sumlist),2))
            CPG_dic[dataline[0]][dataline[1]].append(get_median(sumlist))
    return CPG_dic

filename=sys.argv[1]
# read position of CPG_islands & read position of read_positions (clipping or not) in relation to chr
CPG= parse_file(filename)
f=filename.split("/")[-1].split("-")[0]
chr=f.split("_")[0]
#process across the CPG_islands
CPG_report=get_summary_for_all(CPG)
#print("chr\tnber_CpG_analysed\tCpG_in\tmean_in\tmed_in\tCpG_out\tmean_out\tmed_out\tnber_with_in_prob>out_prob\tnber_with_in_prob<out_prob\tnber_with_in_prob~out_prob\n")
print(chr+"\t"+str(CPG_report[6])+"\t"+str(CPG_report[0])+"\t"+str(round(CPG_report[1]/CPG_report[0],2))+"\t"+str(round(CPG_report[2]/CPG_report[6],2))+"\t"+str(CPG_report[3])+"\t"+str(round(CPG_report[4]/CPG_report[3],2))+"\t"+str(round(CPG_report[5]/CPG_report[6],2))+"\t"+str(CPG_report[8])+"\t"+str(CPG_report[9])+"\t"+str(CPG_report[7]))

