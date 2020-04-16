#!/usr/bin/python

##Copyright (C) 2020 The Institute of Cancer Research (ICR) -- Marco Punta
##
##This Program is free software; you can redistribute it and/or modify it under the
##terms of the GNU General Public License as published by the Free Software Foundation;
##either version 3 of the License, or (at your option) any later version.
##Additional permissions under GNU GPL version 3 section 7:
##This Program is distributed as a service to the research community and is experimental
##in nature and may have hazardous properties.  The Program is distributed WITHOUT ANY WARRANTY,
##express or implied.  In particular all warranties as to SATISFACTORY QUALITY or FITNESS FOR A
##PARTICULAR PURPOSE are excluded. See the GNU General Public License for more details.
##You should have received a copy of the GNU General Public License along with this program;
##if not, see <http://www.gnu.org/licenses>.
##
##You shall not make any use of the name of The Institute of Cancer Research (ICR) in connection
##with this Program in any press or other public announcement without the prior written consent
##of the Institute of Cancer Research.
##Further Restriction: At the request of ICR you will, from time to time provide, free of charge,
##reports on the Program's performance (addressing quality, content and functionality of the Program),
##such reports shall also identify any errors, bugs or shortcomings in the Program as well as the
##Recipient's comments and observations as the ICR may from time to time reasonably request.



#THIS SCRIPT TAKES AS INPUT A FILE WITH A LIST OF PAIRED PRIMARY 
#AND SECONDARY MUTATIONS TOGETHER WITH THE ASSOCIATED WT, PRIMARY
#MUTATION AND SECONDARY MUTATION AMINO ACID SEQUENCES AND GENERATES
#THREE LISTS OF SHORT PEPTIDES:
#LIST OF SHORT PEPTIDES PRESENT IN PRIMARY MUTATION SEQUENCE AND NOT 
#IN THE WT SEQUENCE
#LIST OF SHORT PEPTIDES PRESENT IN SECONDARY MUTATION SEQUENCE AND NOT 
#IN THE WT SEQUENCE
#LIST OF SHORT PEPTIDES PRESENT IN SECONDARY MUTATION SEQUENCE (REVERTANT 
#SEQUENCE) AND NOT IN THE WT SEQUENCE AND NOT IN THE PRIMARY MUTATION SEQUENCE
#THE LENGTH OF SHORT PEPTIDES BEING GENERATED IS DEFINED BY THE USER

import sys
import re

#FUNCTIONS

#TAKES 4 ARGUMENTS: INPUT_FILE,WT_SEQ,PRIMARY_SEQ,REVERTED_SEQ
#RETURNS THREE DICTIONARIES: WT_SEQ,PRIMARY_SEQ,REVERTED_SEQ
def read_sequences(input_seqs_sub,wt_seq_sub,primary_seq_sub,reverted_seq_sub):

    column = {}
    try:
        with open(input_seqs_sub, 'r') as f:
            for i in f:
                ##    case.id    gene    reverted_allele    primary_mutation    secondary_mutation    full_WT_sequence    full_sequence_after_primary_mutation    full_sequence_after_primary_and_secondary_mutation    Query    
                #1    Patient-1-1    BRCA1    c.1360    c.1360    c.1360    MDLSALRV*    MDLSALR*    MDLSA*

                #ASSIGNS TITLE LINE COLUMNS TO DATA
                if re.search("primary_mutation",i):
                    all=i.split("\t")
                    for j in range(0,len(all)):
                        column[all[j]]=int(j)
                    continue

                #ASSIGNS TO EACH MUTATION WT, PRIMARY AND REVERTED SEQUENCES
                i = i.rstrip()
                all=i.split("\t")
                mut=str(all[int(column["gene"])])+"|"+str(all[int(column["primary_mutation"])])+"|"+str(all[int(column["secondary_mutation"])])
                if str(mut) in wt_seq_sub:
                    continue
                wt_seq_sub[mut]=all[int(column["full_WT_sequence"])].split("*")[0]
                primary_seq_sub[mut]=all[int(column["full_sequence_after_primary_mutation"])].split("*")[0]
                reverted_seq_sub[mut]=all[int(column["full_sequence_after_primary_and_secondary_mutation"])].split("*")[0]


    except IOError:
        print input_seqs_sub+" does not exist"


    return(wt_seq_sub,primary_seq_sub,reverted_seq_sub)



#TAKES 1 ARGUMENT: AN AA SEQUENCE 
#RETURNS ONE DICTIONARY WITH ALL PEPTIDES ASSOCIATED TO THE INPUT SEQUENCE
def get_peptides(sub_seq):

    sub_peptides = {}
    aa_seq_list=list(sub_seq)
    stop=len(aa_seq_list)
    aa_set={"A": "1","C": "1","D": "1","E": "1","F": "1","G": "1","H": "1","I": "1","K": "1","L": "1","M": "1","N": "1","P": "1","Q": "1"
,"R": "1","S": "1","T": "1","V": "1","Y": "1","W": "1"}
    
    #GENERATES PEPTIDES FOR ALL K-MER VALUES
    END=KMERS.split(",")
    for j in range(0,len(END)):
        stop_process=0
        for i in range(1,int(stop)+1):
            #STOPS IF AT THE END OF THE SEQUENCE
            if (int(stop)-int(i)+1) < int(END[j]):
                break
            #OTHERWISE, GENERATES PEPTIDE
            pep_seq=[]
            for k in range(0,int(END[j])):
                if str(aa_seq_list[(i-1)+k]) not in aa_set:
                    stop_process=1
                    break
                else:
                    pep_seq.append(aa_seq_list[(i-1)+k])
            
            if int(stop_process) == 1:
                break
            else:
                sub_peptides["".join(pep_seq)]=1


    return(sub_peptides)



        
#MAIN

if len(sys.argv) < 6:
    print "Usage: sys.argv[0] INPUT_sequence_file primary_peps.out rev_peps.out rev_unique_peps.out 8,9,10,11"
    sys.exit()

#INPUT FILE
#SHOULD HAVE A TITLE LINE WITH THE FOLLOWING TAB SEPARATED COLUMNS:
#gene	primary_mutation	secondary_mutation	full_WT_sequence	full_sequence_after_primary_mutation	full_sequence_after_primary_and_secondary_mutation
#full_WT_sequence is the WT AA SEQUENCE OF THE GENE
#full_sequence_after_primary_mutation is the AA SEQUENCE AFTER TAKING INTO ACCOUNT THE PRIMARY MUTATION (e.g. IF PRIMARY 
#MUTATION IS A STOP CODON, SEQUENCE WILL BE TRUNCATED). PRIMARY SEQUENCE HEREAFTER
#full_sequence_after_primary_and_secondary_mutation is the AA SEQUENCE AFTER TAKING INTO ACCOUNT BOTH PRIMARY AND SECONDARY MUTATION.REVERTANT SEQUENCE HEREAFTER
input_seqs=sys.argv[1]
#OUTPUT FILE CONTAINING PEPTIDES FOUND IN THE PRIMARY SEQUENCE AND NOT IN THE WT SEQUENCE
output_primary_peps=sys.argv[2]
#OUTPUT FILE CONTAINING PEPTIDES FOUND IN THE SECONDARY SEQUENCE AND NOT IN THE WT SEQUENCE
output_rev_peps=sys.argv[3]
#OUTPUT FILE CONTAINING PEPTIDES FOUND IN THE SECONDARY SEQUENCE AND NOT IN THE WT SEQUENCE AND NOT IN THE PRIMARY SEQUENCE
output_rev_unique_peps=sys.argv[4]
#LENGTHS OF PEPTIDES BEING PRODUCED (COMMA SEPARATED, e.g. 8,9,10,11)
KMERS=sys.argv[5]

wt_seq = {}
primary_seq = {}
reverted_seq = {}
#read_sequences FUNCTION READS WT, PRIMARY AND REVERTANT SEQUENCE FROM INPUT FILE
(wt_seq,primary_seq,reverted_seq)=read_sequences(input_seqs,wt_seq,primary_seq,reverted_seq)

out_primary_peps=open(output_primary_peps,'w')
out_rev_peps=open(output_rev_peps,'w')
out_rev_unique_peps=open(output_rev_unique_peps,'w')

for mut in wt_seq:
    wt_peptides = {}
    #get_peptides FUNCTION GENERATES ALL PEPTIDES ASSOCIATED TO A GIVEN AA SEQUENCE
    wt_peptides=get_peptides(wt_seq[mut])
    primary_peptides = {}
    primary_peptides=get_peptides(primary_seq[mut])
    reverted_peptides = {}
    reverted_peptides=get_peptides(reverted_seq[mut])
    #WRITES PEPTIDES IN PRIMARY SEQUENCE AND NOT IN WT SEQUENCE
    for pep in primary_peptides:
        if str(pep) not in wt_peptides:
            #CGIFSTAR        N/A     5946delT_S1982Rfs       BRCA2   N/A
            out_primary_peps.write(str(pep)+"\t"+"N/A"+"\t"+str(mut.split("|")[1])+"|"+str(mut.split("|")[2])+"\t"+str(mut.split("|")[0])+"\t"+"N/A"+"\n")
    #WRITES PEPTIDES IN REVERTANT SEQUENCE AND NOT IN WT SEQUENCE
    for pep in reverted_peptides:
        if str(pep) not in wt_peptides:
            out_rev_peps.write(str(pep)+"\t"+"N/A"+"\t"+str(mut.split("|")[1])+"|"+str(mut.split("|")[2])+"\t"+str(mut.split("|")[0])+"\t"+"N/A"+"\n")
    #WRITES PEPTIDES UNIQUE TO REVERTANT SEQUENCE (i.e. NEITHER IN WT NOR IN PRIMARY)
    for pep in reverted_peptides:
        if str(pep) not in wt_peptides and str(pep) not in primary_peptides:
            out_rev_unique_peps.write(str(pep)+"\t"+"N/A"+"\t"+str(mut.split("|")[1])+"|"+str(mut.split("|")[2])+"\t"+str(mut.split("|")[0])+"\t"+"N/A"+"\n")



out_primary_peps.close()
out_rev_peps.close()
out_rev_unique_peps.close()







