#!/usr/bin/env python

from sys import argv
import random
import numpy as np
from scipy.stats import norm
import scipy
import pandas as pd
from itertools import chain 


#alignment.ali has headers with start and end of the domains in the header of fasta separated with ";" like this:
#>EDP05298 pep:known supercontig:v3.1:DS4 ;51;337
#feel free to come up with a better way for passing up that information
#all_phosp file has to contain id, amino acid and number of the phosphoside


#this is to prepare the column names and indexes
def prepare_cols_indx(alignment_file):
	alignment=open(alignment_file,"r").readlines()
	ali_len=len(alignment[1])
	column_names=[]
	indexes=['start1','end1']
	for i in range(0,len(alignment)-1,2):
		protein_name=alignment[i].split(" ")[0].lstrip(">")
		start=alignment[i].split(";")[-2].strip()
		column_names.append(protein_name+"+"+start)

	for i in list(range(0,ali_len-1)):
		indexes.append(i)
	
	return(column_names, indexes)	

#####dataframe with alignment as letters
def letter_ali_dataframe(alignment_file, column_names, indexes):
	alignment=open(alignment_file,"r").readlines()
	ali_len=len(alignment[1])
	alignment_list=[]
	start_row=[]
	end_row=[]
	for i in range(0,len(alignment)-1,2):
		start=int(alignment[i].split(";")[1].strip())
		end=int(alignment[i].split(";")[2].strip())
		start_row.append(start)
		end_row.append(end)

	alignment_list.append(start_row)
	alignment_list.append(end_row)

	for j in range(0,ali_len-1):
		letter_row=[]
		for i in range(1,len(alignment),2):
			letter_row.append(alignment[i][j])
		alignment_list.append(letter_row)
	
	letter_alignment=pd.DataFrame(alignment_list, columns=column_names, index=indexes)
	return(letter_alignment)

############dataframe with alignment as positions
def pos_dataframe(letter_alignment, column_names, indexes):

	position_alignment=pd.DataFrame(index=indexes[2:])

	for protein in column_names:
		pos_seq=[]
		start=letter_alignment.loc['start1',protein]
		seq=letter_alignment[protein].tolist()[2:]
		counter=0
		for aa in seq:
			if aa=="-":
				pos_seq.append(0)
			else:
				pos_seq.append(start+counter)
				counter+=1
		position_alignment[protein]=pos_seq
	return(position_alignment)

###### 0/1 dataframe with mapped phosphorylations onto alignment
def phos_dataframe(domain_name, letter_alignment, position_alignment, column_names, indexes):
	phosp_info=open(domain_name+"_all_phosp","r").readlines()
	
	phosp_alignment=pd.DataFrame(0,index=indexes[2:], columns=column_names)
	for column in column_names:
		protein_id=column.split("+")[0].strip()
		start_no=int(column.split("+")[1].strip())
		for i in range(0,len(phosp_info)):
			l=phosp_info[i].split(",")
			protein_id2=l[1].strip()
			phosp_pos=int(l[2].strip())
			phosp_aa=l[3].strip()
			if protein_id==protein_id2 and start_no<=phosp_pos:
				row_with_pos=position_alignment.index[position_alignment[column]==phosp_pos].tolist()
				if len(row_with_pos)==1:
					phosp_alignment.at[row_with_pos[0],column]+=1
	return(phosp_alignment)

#######rolling window for a list (returns list [2:-2] of original length)
def count_window_for_list(list):                                                 
    window_values=[]                                                             
    for i in range(2,len(list)-2):                                               
        a=int(list[i-2])                                                         
        b=int(list[i-1])                                                         
        c=int(list[i])                                                           
        d=int(list[i+1])                                                         
        e=int(list[i+2])                                                         
        how_many=a+b+c+d+e                                                       
        bg=float(how_many/5.)                                                    
        window_values.append(bg)                                                 
    return(window_values)                                                        
                                                                                 
#background construction
#make_random_histogram(phosp,total_S) of phosphosites across all available S (or T or Y or whatever)
def make_random_histogram(elements, size):
	rlist = [0 for i in range(size)]
	for i in [random.randint(0, size-1) for i in range(elements)]:
		rlist[i] += 1
	return(rlist)

#####################returns one permutated column of rolled window on background
def permutated_dataframe(letter_alignment,phosp_alignment, column_names, indexes):
	permutated_alignment=pd.DataFrame(0, index=indexes[2:], columns=column_names)

	for protein_id in column_names:
		rows_with_S=letter_alignment.index[letter_alignment[protein_id]=="S"].tolist()
		rows_with_T=letter_alignment.index[letter_alignment[protein_id]=="T"].tolist()
		rows_with_Y=letter_alignment.index[letter_alignment[protein_id]=="Y"].tolist()

		phosp_rows=phosp_alignment.index[phosp_alignment[protein_id]>=1].tolist()
		pS=pT=pY=0
		for i in phosp_rows:
			if i in rows_with_S:
				pS+=1
			if i in rows_with_T:
				pT+=1
			if i in rows_with_Y:
				pY+=1
		random_S=make_random_histogram(pS, len(rows_with_S))
		random_T=make_random_histogram(pT, len(rows_with_T))
		random_Y=make_random_histogram(pY, len(rows_with_Y))
		for i,j in zip(random_S, rows_with_S):
			permutated_alignment.at[j,protein_id]=i
		for i,j in zip(random_T, rows_with_T):
			permutated_alignment.at[j,protein_id]=i
		for i,j in zip(random_Y, rows_with_Y):
			permutated_alignment.at[j,protein_id]=i


	sum_of_phosps=permutated_alignment.loc[0:].sum(axis=1).tolist()
	one_of_bg=count_window_for_list(sum_of_phosps)
	return(one_of_bg)

#permutated_alignment=permutated_dataframe(letter_alignment,phosp_alignment, column_names, indexes)
#print(permutated_alignment)
def count_zscore(fg,mean,stdev):
	if stdev==0:
		zscore=0
	else:
		zscore=(fg-mean)/stdev
	
	return(zscore)

def count_pval(fg, mean, stdev):
	fg=float(fg)
	mean=float(mean)
	stdev=float(stdev)
	if fg>mean:
		z=count_zscore(fg, mean, stdev)
		p_val = scipy.stats.norm.sf(abs(z))
	if fg<=mean:
		z=1.
		p_val=1.
	return(p_val)

###############if theres a regulatory file, add the vals for regulatory points
def regulatory_dataframe(domain_name,letter_alignment,position_alignment, column_names, indexes):
	reg_info=open(domain_name+".reg","r").readlines()	

	reg_alignment=pd.DataFrame(0,index=indexes[2:],columns=column_names)
	for column in column_names:
		if "human" in column:
			protein_id=column.split("+")[0].strip().split("#")[0].strip()
			start_no=int(column.split("+")[1].strip())
			for i in range(0,len(reg_info)):
				l=reg_info[i].split(",")
				protein_id2=l[2].strip()
				reg_pos=int(l[1].strip())
				reg_aa=l[0].strip()
				if protein_id==protein_id2 and start_no<=reg_pos:
					row_with_pos=position_alignment.index[position_alignment[column]==reg_pos].tolist()
					if len(row_with_pos)==1:
						reg_alignment.at[row_with_pos[0],column]+=1

	sum_of_regs=reg_alignment.loc[0:].sum(axis=1).tolist()
	return(sum_of_regs[2:-2])

##########active predictions are coming from pfam_scan.pl 
def active_dataframe(domain_name, letter_alignment, position_alignment, column_names, indexes):
	active_info=open(domain_name+".act","r").readlines()
	act_alignment=pd.DataFrame(0,index=indexes[2:],columns=column_names)
	for column in column_names:
		protein_id=column.split("+")[0].strip().split("#")[0].strip()
		start_no=int(column.split("+")[1].strip())
		for i in range(0,len(active_info)):
			protein_id2=active_info[i].split()[0].strip()
			if protein_id==protein_id2:
				active_list=[]
				sites_all=active_info[i].split("predicted_active_site")[1:]
				for k in range(0,len(sites_all)):
					for kk in sites_all[k].split():##this looks extremaly stupid, but the active file is stupid too
						for kkk in kk.split(","):
							active_list.append(int(kkk))
				for a in active_list:
					row_with_pos=position_alignment.index[position_alignment[column]==a].tolist()
					if len(row_with_pos)==1:
						act_alignment.at[row_with_pos[0],column]+=1
	sum_of_acts=act_alignment.loc[0:].sum(axis=1).tolist()
	return(sum_of_acts[2:-2])









###########script calling starts here
#########these are the 3 starting dataframes
alignment_file=argv[1]
domain_name=alignment_file.split(".")[0]
(column_names, indexes)=prepare_cols_indx(alignment_file)
letter_alignment=letter_ali_dataframe(alignment_file, column_names, indexes)
position_alignment=pos_dataframe(letter_alignment, column_names, indexes)
phosp_alignment=phos_dataframe(domain_name, letter_alignment,position_alignment, column_names, indexes)

############foreground (is added below to background_dataframe)
sum_of_phosps=phosp_alignment.loc[0:].sum(axis=1).tolist()                       #
foreground=count_window_for_list(sum_of_phosps)                                  #

####### how many permuts(columns) to create in background dataframe (which are rolled windowed to cound pvals of it)
how_many_permuts=100####################HOW MANY PERMUTS
background_dataframe=pd.DataFrame(index=indexes[4:-2],columns=list(range(0,how_many_permuts)))
for k in range(0,how_many_permuts):
	background_dataframe[k]=permutated_dataframe(letter_alignment,phosp_alignment, column_names, indexes)

##column_names are from 0-100
####
background_dataframe['bg_medians']=background_dataframe.median(axis=1).tolist()
background_dataframe['bg_means']=background_dataframe.mean(axis=1).tolist()
background_dataframe['bg_stdev']=background_dataframe.std(axis=1).tolist()
background_dataframe['foreground']=foreground

#print(background_dataframe)
###############this is to count pvals and add them to the column in background_dataframe
all_pvals=[]
for row in background_dataframe.index:
	fg=background_dataframe.at[row,'foreground']
	mean=background_dataframe.at[row,'bg_means']
	st_dev=background_dataframe.at[row,'bg_stdev']
	p_val=count_pval(fg, mean, st_dev)
	all_pvals.append(p_val)

background_dataframe['pvals']=all_pvals
background_dataframe['pvals']=background_dataframe['pvals']+0.00000000000000001 ###############in case some values can be 0
background_dataframe['minus_logs']=-np.log10(background_dataframe.pvals)
background_dataframe['regs']=regulatory_dataframe(domain_name, letter_alignment, position_alignment, column_names, indexes)
background_dataframe['active']=active_dataframe(domain_name, letter_alignment, position_alignment, column_names, indexes)
header=["foreground", "bg_means", "bg_stdev", "pvals","minus_logs", "regs","active"]

#############write down datagrame
background_dataframe.to_csv('%s_pval_full_ali'%domain_name,columns=header)
