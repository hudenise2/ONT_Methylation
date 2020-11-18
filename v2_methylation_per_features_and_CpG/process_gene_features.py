#!/usr/bin/python3

import argparse, re
import sys, statistics as st

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
           print(line)
           print(prefix)
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
                    feature_dic['5mRNA'][mRNA_name][1]=start
                 else:
                    feature_dic['3mRNA'][mRNA_name][1]=start
    return feature_dic

def extract_CPG_per_feature(CPG, features):
    '''
    Function to associate CpG to features
    CPG {start : [end, name], ...} such as {1:[912, 'CpG1'], 1429:[2152, 'CpG2'], ...}
    features {'5gene': {'ENSACLG00000000003': [78636, 79235], 'ENSACLG00000000004': [85876, 85876],...}, '3gene': {'ENSACLG00000000002': [19687, 19687]}, '5mRNA': {'ENSACLT00000000003': [78636, 80433], 'ENSACLT00000000004': [79235, 83444],...}, '3mRNA': {'ENSACLT00000000002': [19687, 21573]}, 'CDS': {'ENSACLP00000000002': [19687, 21631], 'ENSACLP00000000003': [79235, 80801],...}
    '''
    res_dic={}
    end_CPG={y[0]:y[1] for x, y in CPG.items()}
    for feature in features:
        if feature not in res_dic:
           res_dic[feature]={}
        for entry in features[feature]:
           start = min(features[feature][entry])
           end = max(features[feature][entry])
           entry_length=end - start
           if entry_length!=0:
              res_dic[feature][entry]=[entry_length]
              CpG_island_list=[CPG[x][1] for x in CPG.keys() if int(x) >= start and int(x) <= end]
              end_CpG_island_list=[y[1] for x,y in CPG.items() if int(x) >= start and int(x) <= end]
              final_CpG_island_list=list(set(CpG_island_list+end_CpG_island_list))
              res_dic[feature][entry].append(final_CpG_island_list)
           elif feature == 'between': 
              '''
              This is to deal with the case where the between feature goes to the end of the chromosome.
              The gene_feature files I used  do not contain the position of the end of the chromosome... I had to extract them from the gff file
              at ~/rds/rds-rd109-durbin-group/ref/fish/Astatotilapia_calliptera/fAstCal1.2/gff3/Astatotilapia_calliptera.fAstCal1.2.99.primary_assembly.*.gff3.gz
              '''
              end_chrom={'chr1':41162407,'chr2':38215696,'chr3':52512415,'chr4':33433079,'chr5':38678279,'chr6':41428020,'chr7':68397778,'chr8':25827708,'chr9':35913139,'chr10':34074121,'chr11':36676067,'chr12':38669361,'chr13':33145951,'chr14':39985354,'chr15':40222765,'chr16':35801853,'chr17':39448915,'chr18':35852588,'chr19':30863130,'chr20':31467755,'chr22':35011464,'chr23':43921804}
              res_dic[feature][entry]=[end_chrom[prefix]-start]
              CpG_island_list=[CPG[x][1] for x in CPG.keys() if int(x) >= start]
              end_CpG_island_list=[y[1] for x,y in CPG.items() if int(x) >= start]
              final_CpG_island_list=list(set(CpG_island_list+end_CpG_island_list))
              res_dic[feature][entry].append(final_CpG_island_list)
    return res_dic
        
def process_prob(CPG_dic, prob_dic):
    high_prob={}   #will contain the mean probability/CpG and length for the entries per feature having this probability/CpG > 192
    prob_feature_dic={}  #will contain the mean probability per CpG and length for each entry per feature  
    all_results={}    #will contain list of length, probabiliy of Cm and number of CpG per feature
    for feature in CPG_dic:
       CpG_all_list=[]   #will contain the number of CpG for each entry per feature
       prob_all_list=[]  #will contain the sum of the Cm probability for each entry per feature
       feature_length_list=[]  #will contain the length each entry per feature
       if feature not in high_prob:
          high_prob[feature]={}
       if feature not in prob_feature_dic:
          prob_feature_dic[feature]={}
       if feature not in all_results:
          all_results[feature]=[]
       for entry in CPG_dic[feature]:
          feature_length_list.append(CPG_dic[feature][entry][0])
          #list of CpG islands for the entry in feature
          CpG_list = CPG_dic[feature][entry][1]
          #list of Cm probability for the entry in feature 
          prob_list = [int(prob_dic[x]) for x in CpG_list if x in prob_dic]
          #calculate the probability across the probability in feature
          prob_sum=sum(prob_list)
          prob_all_list.append(prob_sum)
          #number of CPG per feature
          count = len(CpG_list)
          CpG_all_list.append(count)
          if count > 0:
             #get the features with high probabilities (defined as having probability/number of CpG islands > 192 i.e. each CpG islands have on average prob > 192)
             if prob_sum/count > 192:
                high_prob[feature][entry]=[round(float(prob_sum/count),1), CPG_dic[feature][entry][0]] 
             prob_feature_dic[feature][entry]=[round(float(prob_sum/count),1), CPG_dic[feature][entry][0]]
          else: prob_feature_dic[feature][entry]=[0, CPG_dic[feature][entry][0]]
       all_results[feature]=[feature_length_list, prob_all_list, CpG_all_list]
       print("FEATURE: "+feature)
       print("   number of feature: "+str(len(CPG_dic[feature])))
       if len(CPG_dic[feature])>0:
          print("   mean length of feature (min, max, median): "+str(round(st.mean(feature_length_list),1))+" ("+str(min(feature_length_list))+", "+str(max(feature_length_list))+", "+str(round(st.median(feature_length_list),1))+")")
          print("   features with no CPG: "+str(len(CPG_dic[feature]) - len({x:y for x, y in prob_feature_dic[feature].items() if y[0] > 0})))
          print("   features with CPG: "+str(len({x:y for x, y in prob_feature_dic[feature].items() if y[0] > 0})))
          print("   mean number of CpG by feature (min, max, median): "+str(round(st.mean(CpG_all_list),1))+" ("+str(min(CpG_all_list))+", "+str(max(CpG_all_list))+", "+str(round(st.median(CpG_all_list),1))+")")
          if len(prob_feature_dic[feature]) >0 and sum(prob_all_list) > 0:
             print("   mean probability per feature (min, max, median): "+str(round(st.mean(prob_all_list),1))+" ("+str(min(prob_all_list))+", "+str(max(prob_all_list))+", "+str(round(st.median(prob_all_list),1))+")")
          #list of probability per feature and CpG
          list_prob_per_CpG = [y[0] for x, y in prob_feature_dic[feature].items() if y[0] > 0]
          print("   mean probability per feature/CpG (min, max, median): "+str(round(st.mean(list_prob_per_CpG),1)) +" ("+str(min(list_prob_per_CpG))+", "+str(max(list_prob_per_CpG))+", "+str(round(st.median(list_prob_per_CpG),1))+")")
          #list of prob per feature and length
          prob_per_length = sum(prob_all_list)/float(sum(feature_length_list)/100)
          print("   mean probability per 100 nt: "+str(round(prob_per_length,1)))
          print("   number of feature with high probability: "+str(len(high_prob[feature]))+"\n")
    return high_prob, all_results, prob_feature_dic
           
def parse_output_CPG_Threshold(name):
    sum_dic={}
    for line in open("Threshold"+name+"_output_parse_dic","r"):
       line=line.rstrip()
       dataline=line.split()
       sum=0
       if dataline[1]=="ON":
          for val in dataline[2][1:-1].split(","):
              if len(val) > 0 : sum+=int(val)
          if dataline[0] not in sum_dic:
              sum_dic[dataline[0]]=[sum, len(dataline[2][1:-1].split(","))]
          else:
              sum_dic[dataline[0]][0]+=sum
              sum_dic[dataline[0]][1]+=len(dataline[2][1:-1].split(","))
    prob_dic={k:round(v[0]/v[1],1) for k,v in sum_dic.items()}
    return prob_dic
    

def parse_CPGfile(name):
    CPG_dic={}
    sumprob=0
    medianval=0
    for line in open("CpG_"+name, "r"):
        line=line.rstrip()
        dataline=line.split()
        CPG_dic[dataline[0]] = [int(dataline[1]), dataline[2]]
    return CPG_dic

#gene feature file
filename=sys.argv[1]

#get prefix for other files to process
prefix=filename.split("/")[-1].split("_")[0].replace("omosome","")

# read name of CPG_islands & position. The returned dictionary will be in this format: {CpG2865': [41143415, 41144464], ...
CPG_dic= parse_CPGfile(prefix)
# open the gene feature file. The returned format will be: {1: ['b', 'start'], 19687: ['c', 'ENSACLP00000000002'], 21631: ['5', 'ENSACLP00000000002'], ...
features_dic=parse_genefile(filename)
# get the average Cm prob for the CPG_islands. Format:   {'CpG2834': 9.0, 'CpG2835': 16.7,
prob_dic=parse_output_CPG_Threshold(prefix)
#find CpG islands associated per feature
res_dic=extract_CPG_per_feature(CPG_dic, features_dic)
#associate average Cm probability with features
High_prob, list_results, prob_CpG_dic=process_prob(res_dic, prob_dic)

outputH=open("gene_features/High_prob/CPG_feature_"+prefix+"_results", "w")
#get number and prob per feature
feature_order=['5gene', '5mRNA', 'CDS', '3mRNA', '3gene', 'between']
feature_eq={'5gene':'PROMOTER', '3gene':'END GENE', 'CDS': 'CODING REGIONS', '5mRNA': "5' UTR", '3mRNA': "3' UTR", 'between': "BETWEEN GENES"}
for feat in feature_order:
    outputH.write(feature_eq[feat]+"\n")
    outputH.write(str(High_prob[feat])+"\n")
outputH.close()

output=open("gene_features/CPG_feature_"+prefix+"_results", "w")
output.write("feature\ttotal_number\tfeatures_with_noCPG\tfeatures_with_CPG\tmean_number_CPG/feature\tmean_length/feature\tmean_prob/feature\tmedian_prob/feature\tmax_prob_(nber)\tmin_prob_(nber)\tnumber_with_prob>192\n")
for feat in list_results:
   line=[feature_eq[feat]]
   total=len(list_results[feat][0])
   line.append(str(total))
   if total > 0:
      line.append(str(list_results[feat][2].count(0))+"("+str(round(list_results[feat][2].count(0)*100/total,1))+"%)")
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
