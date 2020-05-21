#############################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies
###
### CITATION:
###
### PROCESS: To obtain the number of clones per sample
###
### Author: Silvia Pineda
### Date: January, 2017
############################################################################################
###import pandas library to work with data frames
import pandas as pd
import math
from functools import partial
from multiprocessing import Pool
################################
######### Functions ############
################################

def MatchNucleotides( a,b ):
    misMatches = 0
    max_mismatches = int(len(a) * 0.10)
    for i in range(0, len(a)):
        if(a[i]!=b[i]):
            misMatches = misMatches + 1
            if misMatches > max_mismatches:
                return False
    return True

def NucleotideBelongsClone(nucleotide,clone):
    for i in range(0,len(clone)):
        matches=MatchNucleotides(nucleotide,clone[i])
        if(matches):
            return True
    return False

def AddNucleotide(nucleotide,clones):
    result_clones=[]
    nucleotide_clone = [nucleotide]
    for i in range(0,len(clones)):
        currentClone = clones[i]
        if(NucleotideBelongsClone(nucleotide,currentClone)):
            nucleotide_clone = nucleotide_clone + currentClone
        else:
            result_clones.append(currentClone)
    result_clones.append(nucleotide_clone)
    return result_clones

def ProcessNucleotides(nucleotides_Array):
    Listclones=[]
    for i in range(0, len(nucleotides_Array)):
        Listclones = AddNucleotide(nucleotides_Array[i], Listclones)
    return AddIndexToClones(Listclones)
    #return Listclones
    
def AddIndexToClones(clones_Array):
    df_clones = pd.DataFrame([])
    for i in range(0, len(clones_Array)):
        cloneNumber = i + 1
        clone = clones_Array[i]
        for j in range(0, len(clone)):
            df_clones = df_clones.append({'clone':clone[j],
                                          'number':cloneNumber},ignore_index=True)
    return df_clones

def ProcessGroup(nucleotides_Dataframe, unique_VJCDR3_Series, i):
    print (i)
    selected_VJCDR3_Dataframe = nucleotides_Dataframe[nucleotides_Dataframe['V_J_lenghCDR3'] == unique_VJCDR3_Series[i]]
    nucleotides_Array = list(selected_VJCDR3_Dataframe['nSeqCDR3'])
    ##Obtain the number of clones inferred
    Clones_Infered_indexed = ProcessNucleotides(nucleotides_Array)
    selected_VJCDR3_Dataframe_sorted = selected_VJCDR3_Dataframe.sort_values(['nSeqCDR3'])
    Clones_Infered_indexed_sorted = Clones_Infered_indexed.sort_values(['clone'])
    selected_VJCDR3_Dataframe_sorted['CloneId'] = list(Clones_Infered_indexed_sorted['number'])
    return selected_VJCDR3_Dataframe_sorted

def ProcessSample(nucleotides_Dataframe):
    result = pd.DataFrame([])
    unique_VJCDR3_Series = nucleotides_Dataframe['V_J_lenghCDR3'].unique()
    ##parallelization for each unique VJCDR3
    pool = Pool(30)
    partialFunction = partial(ProcessGroup, nucleotides_Dataframe, unique_VJCDR3_Series)
    result = pool.map(partialFunction, range(0,len(unique_VJCDR3_Series)))
    return result

def main():
    ########################
    ###### Main program ####
    ########################
    print("Start")
    nucleotides_Dataframe = pd.read_csv("~/TCGA-Immune/Data/Validation_Normal_pancreas/data_for_cloneInfered_TCR_PancreasValidation.txt",sep="\t")
    
    result_ClonesInfered = pd.DataFrame([])
    result = ProcessSample(nucleotides_Dataframe)
    result_ClonesInfered = result_ClonesInfered.append(result)

    ###Result
    result_ClonesInfered.to_csv('~/TCGA-Immune/Data/Validation_Normal_pancreas/ClonesInfered_TCR_ValidationNormal.csv')
    print("End")

main()
