#!/usr/bin/python3

import argparse, re, sys, logging 
from datetime import datetime
import statistics as st
import pysam
import random
from Bio import SeqIO


logging.basicConfig(filename='proc_e',  level=logging.DEBUG)

def parse_genefile(file):
    '''
    function to parse the file gene_features/chromosome<X>_features.
    The format of this file is <feature_type> <feature_name> <strand> <start_po> <end_pos>
    The output format is a dictionary with feature as keys and nested dictionary with feature_id as keys and list of start and end position for each feature:
      {'5gene' : {'gene_id1':[start, end], 'gene_id2':[start, end]... }, '5mRNA':{'mRNA1':[start, end], 'mRNA2':[start, end],...}, ...}
    model used:
            |--------------------------gene-----------------------------|
      model .   |-------------------mRNA----------------------------|   .
            .   .   |-5UTR-|  |--------CDS-----------|  |-3UTR-|    .   .
            .   .             .                      .              .   .
             ___ ____________ _______________________ ______________ ____...
     feature 5gene   5mRNA             CDS                3mRNA      3gene b  

    '''
    feature_dic={}
    feature_dic['5gene']={}
    feature_dic['3gene']={}
    feature_dic['5mRNA']={}
    feature_dic['3mRNA']={}
    feature_dic['CDS']={}
    between_count=0
    feature_dic["between"]={}
    feature_dic["between"][between_count]=[1, 1]
    for line in open(file, "r"):
        line=line.rstrip()
        dataline=line.split()
        strand=dataline[-3]
        start=int(dataline[-2])
        end=int(dataline[-1])
        name=dataline[1].split(";")[0].split(":")[1]
        cds_flag=0
        #only considering 3 types of features
        if dataline[0] in ['gene', 'mRNA', 'CDS']:
           if dataline[0] == 'gene':
              if strand=="+":
                 feature_dic['5gene'][name]=[start, end]
              else:
                 feature_dic['3gene'][name]=[start, end]
              between_count+=1
              feature_dic["between"][between_count]=[end, end]
              feature_dic["between"][between_count - 1][0]=start   
           if dataline[0]== 'mRNA':
              gene_name=dataline[1].split(";")[1].split(":")[1]
              if strand=="+":
                 feature_dic['5gene'][gene_name][1]=start 
                 feature_dic['5mRNA'][name]=[start, end]
              else:
                 feature_dic['3gene'][gene_name][1]=start 
                 feature_dic['3mRNA'][name]=[start, end]
           if dataline[0]== 'CDS':
              mRNA_name=dataline[1].split(";")[1].split(":")[1]
              if name in feature_dic['CDS']:
                 feature_dic['CDS'][name][1]=start
              else:
                 feature_dic['CDS'][name]=[start, end]
                 if strand=="+":
                    feature_dic['5mRNA'][mRNA_name][1]=end
                 else:
                    feature_dic['3mRNA'][mRNA_name][1]=start
    return feature_dic

def get_reads_per_feature(feature_dic, mod_dic, pos_dic, CpG_dic):
   '''
   this script will return a dictionary  ex {'5gene':{'name_read1:{'noCpG':[list of CG probability],'inCpG': [list of CG probability]}, {'name_read2:{'noCpG':[list of CG probability],'inCpG': [list of CG probability]}, ...}}, {'5mRNA':{....}, ...}}
   '''
   prob_res={}
   #the line below is only used for testing
   feature_dic={'5gene': {'ENSACLG00000020698': [[32873740, 32874604], []], 'ENSACLG00000008971': [[23856159, 23856667], ['CpG1541']]}, '3gene': {'ENSACLG00000006950': [[4463837, 4477631], ['CpG404','CpG403']]}, '5mRNA': {'ENSACLT00000000003': [[78636, 80433], []], 'ENSACLT00000000010': [[155986, 162214], ['CpG18', 'CpG19']]}, '3mRNA': {'ENSACLT00000015888': [[6373039, 6375498], []]}, 'CDS': {'ENSACLP00000000002': [[19687, 21573], ['CpG5']], 'ENSACLP00000000003': [[79235, 80433], ['CpG10']]}, 'between': {1:[[21631, 78636], ['CpG8', 'CpG9', 'CpG6', 'CpG7']], 2:[[83577, 85876], ['CpG11']]}}
   bamfile = pysam.AlignmentFile("/home/had38/rds/rds-rd109-durbin-group/projects/cichlid/ONT_Basecall/Hubert_OtoArg_fin_L12_2020-03-16/Align_PASS/fastq_mapping_sorted.bam", "rb")
   for feat in feature_dic:
      if feat not in prob_res:
        prob_res[feat]={}
      for entry in feature_dic[feat]:
         prob_res[feat][entry]={}
         start=feature_dic[feat][entry][0][0]
         end=feature_dic[feat][entry][0][1]
         CpG_list=feature_dic[feat][entry][1]
         #get coordinate of CpG islands
         CpG_coord_list=[[CpG_dic[x][0], CpG_dic[x][1]] for x in CpG_list]
         logging.info("\n"+str(entry))
         logging.info(datetime.now())
         iter = bamfile.fetch("chr1", start, end)
         logging.info(datetime.now())
         prob_res[feat][entry]['length']=0
         for read in iter:
            seq_name=str(read.query_name)
            prob_res[feat][entry][seq_name]={}
            seq_align_start= int(read.query_alignment_start) #read coord
            seq_align_end = int(read.query_alignment_end) #read coord
            chr_align_start= int(read.reference_start) #corresponding chr coord
            chr_align_end = int(read.reference_end) #corresponding chr coord
            read_length=chr_align_end -  chr_align_start
            align_flag=int(read.flag)
            seq = str(read.query_alignment_sequence) #just the portion of the sequence which is aligned
            #get position of CG dinucleotide in seq
            #if align_flag==16 : CG_pos=[(x.start()+1) for x in re.finditer("GC", seq)] # read coord. Reverse complement aligned
            #elif align_flag==0: CG_pos=[x.start() for x in re.finditer("CG", seq)] # read coord. 
            #CG_pos=pos_dic[seq_name]
            #get probability for these positions
            #prob=[int(mod_dic[seq_name][x]) for x in CG_pos]
            #prob=[int(random.uniform(0, 255)) for x in CG_pos]
            if seq_name in pos_dic:
               all_CG_pos=pos_dic[seq_name]
               prob=[mod_dic[seq_name][x] for x in range(len(mod_dic[seq_name])) if pos_dic[seq_name][x] >= seq_align_start and pos_dic[seq_name][x] <= seq_align_end]
               CG_pos=[x for x in all_CG_pos if x >= seq_align_start and x <= seq_align_end]
               #extract probability correspond to CG in CpG islands and then the remaining sequences of the features
               if len(CpG_list) > 0:
                 CpG_allCG_index=[]
                 for CpG_index in range(len(CpG_coord_list)):
                    start_CpG=CpG_coord_list[CpG_index][0] #chr coord
                    end_CpG=CpG_coord_list[CpG_index][1]  #chr coord
                    length_CpG=end_CpG - start_CpG
                    read_start_CpG=start_CpG -chr_align_start  #transform in read coord
                    read_end_CpG= end_CpG - chr_align_start   #transform in read coord
                    #get the index of the CG in the CpG island
                    CpG_CG_index=[x for x in range(len(prob)) if CG_pos[x] >= read_start_CpG and CG_pos[x] <= read_end_CpG] 
                    CpG_allCG_index+=CpG_CG_index
                    #get the probability corresponding to the CpG index
                    prob_CpG_list=[prob[x] for x in CpG_CG_index]
                    if 'inCPG' in prob_res[feat][entry][seq_name]:
                       prob_res[feat][entry][seq_name]['inCpG'][0][0]+=prob_CpG_list
                       prob_res[feat][entry][seq_name]['inCpG'][0][1]+=length_CpG
                    else:
                       prob_res[feat][entry][seq_name]['inCpG']=[prob_CpG_list]
                       prob_res[feat][entry][seq_name]['inCpG'].append(length_CpG)                  
                 noCpG_CG_index=[x for x in range(len(CG_pos)) if x not in CpG_allCG_index] 
                 prob_res[feat][entry][seq_name]['noCpG']=[[prob[x] for x in noCpG_CG_index]]
                 prob_res[feat][entry][seq_name]['noCpG'].append(read_length - prob_res[feat][entry][seq_name]['inCpG'][1])
               else: 
                 prob_res[feat][entry][seq_name]['noCpG']=[prob]
                 prob_res[feat][entry][seq_name]['noCpG'].append(read_length)
                 prob_res[feat][entry][seq_name]['inCpG']=[[],0]
            else:
               prob_res[feat][entry][seq_name]['inCpG']=[[],0]
               prob_res[feat][entry][seq_name]['noCpG']=[[],0]
         prob_res[feat][entry]['length']=end - start
   return prob_res

def parse_mod():
    c=0
    mod_dic={}
    gc_pos_dic={}
    logging.info(datetime.now())
    with open("/home/had38/rds/rds-rd109-durbin-group/projects/cichlid/ONT_Basecall/Hubert_OtoArg_fin_L12_2020-03-16/modif_tables/Preprocessed_Table_mod2", "r") as probfile:
       for line in probfile:
           line=line.rstrip()
           if line.startswith(">"):
              seqname=line[1:] 
           elif len(line) > 0:
              mod_dic[seqname]=[int(x) for x in line.split(",")]
    logging.info(datetime.now())
    with open("/home/had38/rds/rds-rd109-durbin-group/projects/cichlid/ONT_Basecall/Hubert_OtoArg_fin_L12_2020-03-16/modif_tables/Preprocessed_Table_pos") as posfile:
       for line in posfile:
           line=line.rstrip()
           if line.startswith(">"):
              seqname=line[1:]         
           elif len(line) > 2:
              gc_pos_dic[seqname]=[int(x) for x in line.split(",")]
    logging.info(datetime.now())
    return mod_dic, gc_pos_dic

def extract_CpG_per_feature(features, CpG):
    '''
    Function to associate CpG to features
    CpG {name : [start, end], ...} such as {'CpG1':[1,912], 'CPG2':[1421, 2152], ...}

    features {'5gene': {'ENSACLG00000000003': [78636, 79235], 'ENSACLG00000000004': [85876, 85876],...}, '3gene': {'ENSACLG00000000002': [19687, 19687]}, '5mRNA': {'ENSACLT00000000003': [78636, 80433], 'ENSACLT00000000004': [79235, 83444],...}, '3mRNA': {'ENSACLT00000000002': [19687, 21573]}, 'CDS': {'ENSACLP00000000002': [19687, 21631], 'ENSACLP00000000003': [79235, 80801],...}
    '''
    res_dic={}
    start_CpG={y[0]:x for x, y in CpG.items()}
    end_CpG={y[1]:x for x, y in CpG.items()}
    for feature in features:
        if feature not in res_dic:
           res_dic[feature]={}
        for entry in features[feature]:
           start = min(features[feature][entry])
           end = max(features[feature][entry])
           entry_length=end - start
           if entry_length!=0:
              res_dic[feature][entry]=[[start, end]]
              CpG_island_list=[start_CpG[x] for x in start_CpG.keys() if int(x) >= start and int(x) <= end]
              end_CpG_island_list=[end_CpG[x] for x in end_CpG.keys() if int(x) >= start and int(x) <= end]
              final_CpG_island_list=list(set(CpG_island_list+end_CpG_island_list))
              res_dic[feature][entry].append(final_CpG_island_list)
           elif feature == 'between': 
              '''
              This is to deal with the case where the between feature goes to the end of the chromosome.
              The gene_feature files I used  do not contain the position of the end of the chromosome... I had to extract them from the gff file
              at ~/rds/rds-rd109-durbin-group/ref/fish/Astatotilapia_calliptera/fAstCal1.2/gff3/Astatotilapia_calliptera.fAstCal1.2.99.primary_assembly.*.gff3.gz
              '''
              end_chrom={'chr1':41162407,'chr2':38215696,'chr3':52512415,'chr4':33433079,'chr5':38678279,'chr6':41428020,'chr7':68397778,'chr8':25827708,'chr9':35913139,'chr10':34074121,'chr11':36676067,'chr12':38669361,'chr13':33145951,'chr14':39985354,'chr15':40222765,'chr16':35801853,'chr17':39448915,'chr18':35852588,'chr19':30863130,'chr20':31467755,'chr22':35011464,'chr23':43921804}
              res_dic[feature][entry]=[[start, end_chrom[prefix]]]
              CpG_island_list=[start_CpG[x] for x in start_CpG.keys() if int(x) >= start]
              end_CpG_island_list=[end_CpG[x] for x in end_CpG.keys() if int(x) >= start]
              final_CpG_island_list=list(set(CpG_island_list+end_CpG_island_list))
              res_dic[feature][entry].append(final_CpG_island_list)
    return res_dic
        
def summary_prob(results_dic):
    '''
    This script will return 2 dictionaries:
    feature_entry summarize the data from each entry for each feature. The lists for 'inCpG' and 'noCpG' contains the sum of probability, the number of CG and the length of the feature.
     ex: {'5mRNA': {'ENSACLT00000000003': {'inCpG': [0, 0, 0], 'noCpG': [4872705, 38398, 27547]}, 'ENSACLT00000000010': {'inCpG': [93264, 688, 3250], 'noCpG': [2929664, 23017, 723986]}}, '3mRNA': {...}, ...} 
    summary_entry summarize the data per feature. ex: {'5mRNA': {'inCpG': [86020, 688, 7429], 'noCpG': [7796966, 61415, 612003]}, '3mRNA': {'inCpG': [0, 0, 0], 'noCpG': [3108491, 24510, 254628]}, 'CDS': {...}, ...}
    '''
    summary_entry={}
    feature_entry={}
    for feat in results_dic:
       feature_entry[feat]={}
       summary_entry[feat]={}
       for entry in results_dic[feat]:
          feature_entry[feat][entry]={}
          count_CG_noCpG=0
          sum_prob_noCpG =0
          len_noCpG=0
          count_CG_inCpG=0
          sum_prob_inCpG=0
          len_inCpG=0
          for read in results_dic[feat][entry]:
             if read == 'length':
                summary_entry[feat]['length']=results_dic[feat][entry]['length']
             else:
                if len(results_dic[feat][entry][read]['noCpG'][0]) > 0:
                   count_CG_noCpG += len(results_dic[feat][entry][read]['noCpG'][0])
                   sum_prob_noCpG += sum(results_dic[feat][entry][read]['noCpG'][0])
                   len_noCpG+= results_dic[feat][entry][read]['noCpG'][1]
                if len(results_dic[feat][entry][read]['inCpG'][0]) > 0:
                   count_CG_inCpG += len(results_dic[feat][entry][read]['inCpG'][0])
                   sum_prob_inCpG += sum(results_dic[feat][entry][read]['inCpG'][0])
                   len_inCpG+= results_dic[feat][entry][read]['inCpG'][1]
          feature_entry[feat][entry]['noCpG']=[sum_prob_noCpG/len(results_dic[feat][entry]), count_CG_noCpG /len(results_dic[feat][entry]), len_noCpG/len(results_dic[feat][entry])]
          feature_entry[feat][entry]['inCpG']=[sum_prob_inCpG/len(results_dic[feat][entry]), count_CG_inCpG/len(results_dic[feat][entry]), len_inCpG/len(results_dic[feat][entry])]
       #return average probability per feature
       summary_entry[feat]['noCpG']=[sum([feature_entry[feat][x]['noCpG'][0] for x in feature_entry[feat]]) / len(feature_entry[feat]), sum([feature_entry[feat][x]['noCpG'][1] for x in feature_entry[feat]])/len(feature_entry[feat]), sum([feature_entry[feat][x]['noCpG'][2] for x in feature_entry[feat]])/len(feature_entry[feat])]
       summary_entry[feat]['inCpG']=[sum([feature_entry[feat][x]['inCpG'][0] for x in feature_entry[feat]]) / len(feature_entry[feat]), sum([feature_entry[feat][x]['inCpG'][1] for x in feature_entry[feat]])/len(feature_entry[feat]), sum([feature_entry[feat][x]['inCpG'][2] for x in feature_entry[feat]])/len(feature_entry[feat])]
    return summary_entry, feature_entry

def parse_CpGfile(name):
    CpG_dic={}
    sumprob=0
    medianval=0
    for line in open("/home/had38/rds/rds-rd109-durbin-group/projects/cichlid/ONT_Basecall/Hubert_OtoArg_fin_L12_2020-03-16/CpG_islands_per_chr/CpG_"+name, "r"):
        line=line.rstrip()
        dataline=line.split()
        CpG_dic[dataline[2]] = [int(dataline[0]), int(dataline[1])]
    return CpG_dic

#gene feature file
filename=sys.argv[1]

#get prefix for other files to process
prefix=filename.split("/")[-1].split("_")[0].replace("omosome","")

# read name of CpG_islands & position. The returned dictionary will be in this format: {CpG2865': [41143415, 41144464], ...
CpG_dic= parse_CpGfile(prefix)
# open the gene feature file. The returned format will be: {'5gene': {'ENSACLG00000000003': [78636, 79235],...}, '3gene': {'ENSACLG00000000002': [19687, 19687]}, '5mRNA': {'ENSACLT00000000003': [78636, 80433],...}, '3mRNA': {'ENSACLT00000000002': [19687, 21573]}, 'CDS': {'ENSACLP00000000002': [19687, 21631],...}
features_dic=parse_genefile(filename)
# associate the CpG to the features. Will return a nested dictionary: {features:{feature_name:[[feature_start, feature_end], [list of CpG islands]],...},...}
CpG_features_dic=extract_CpG_per_feature(features_dic, CpG_dic)
# get the probability and position for each CG position of the reads
mod_dic, pos_dic=parse_mod()
# extract the reads corresponding to each gene features from the sorted bam file and associate the probability for each CG in the feature (in and out of the CpG islands)
raw_results= get_reads_per_feature(CpG_features_dic, mod_dic, pos_dic, CpG_dic)
summary_results, feature_results=summary_prob(raw_results)

#formating and ordering
feature_order=['5gene', '5mRNA', 'CDS', '3mRNA', '3gene', 'between']
feature_eq={'5gene':'PROMOTER', '3gene':'END GENE', 'CDS': 'CODING REGIONS', '5mRNA': "5' UTR", '3mRNA': "3' UTR", 'between': "BETWEEN GENES"}

#outputH=open("gene_features/High_prob/CpG_feature_"+prefix+"_results", "w")
for feat in feature_order:
    #outputH.write(feature_eq[feat]+"\n")
    print(feature_eq[feat])
    print(feature_results[feat])
    print([[x, feature_results[feat][x]['inCpG'][0]/feature_results[feat][x]['inCpG'][1]] for x in feature_results[feat] if int(feature_results[feat][x]['inCpG'][1]) > 0])
    print([[x, feature_results[feat][x]['noCpG'][0]/feature_results[feat][x]['noCpG'][1]] for x in feature_results[feat] if int(feature_results[feat][x]['noCpG'][1]) > 0])
    #list_CpG_high =[x for x in feature_results[feat] if feature_results[feat][x]['inCpG'][1] > 0 and feature_results[feat][x]['inCpG'][0]/feature_results[feat][x]['inCpG'][1] > 191]
    #print(list_CpG_high)
    #outputH.write(str(High_prob[feat])+"\n")
#outputH.close()
'''
#{'5mRNA': {'inCpG': [86020, 688], 'noCpG': [7796966, 61415]}, '3mRNA': {'inCpG': [0, 0], 'noCpG': [3108491, 24510]}, 'CDS':
output=open("gene_features/CpG_feature_"+prefix+"_results", "w")
output.write("feature\ttotal_number\tfeatures_with_noCpG\tfeatures_with_CpG\tmean_number_CpG/feature\tmean_length/feature\tmean_prob/feature\tmedian_prob/feature\tmax_prob_(nber)\tmin_prob_(nber)\tnumber_with_prob>192\n")
for feat in feature_order:
   line=[feature_eq[feat]]
   total=len(summary_results[feat])
   line.append(str(total))
   if total > 0:
      
      line.append(str(summary_results[feat]['noCPG'].count(0))+"("+str(round(list_results[feat][2].count(0)*100/total,1))+"%)")i
      with_CpG=total - list_results[feat][2].count(0)
      line.append(str(with_CpG)+"("+str(round(with_CpG*100/total,1))+"%)")
      line.append(str(round(st.mean(list_results[feat][2]),1)))
      line.append(str(round(st.mean(list_results[feat][0]),1)))
      line.append(str(round(st.mean(list_results[feat][1]),1)))
      line.append(str(round(st.median(list_results[feat][1]),1)))
      line.append(str(max(list_results[feat][1])) + "("+str(list_results[feat][1].count(max(list_results[feat][1])))+")")
      line.append(str(min(list_results[feat][1])) + "("+str(list_results[feat][1].count(min(list_results[feat][1])))+")")
      line.append(str(len(High_prob[feat])))
   output.write("\t".join(line)+"\n")
'''
