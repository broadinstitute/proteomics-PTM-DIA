import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os
import matplotlib.pyplot as plt

#For TV, but can be used to read in a file with spiked peptides as a column
lights = list(pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Lights.tsv', delimiter= '\t')['Modified'])       #Modified sequence annotates phosphopeptides as they appear on Spectronaut
heavies = list(pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Heavies.tsv', delimiter= '\t')['LibraryHeavies'][0:233])

#USER INPUT: Assign variables here#
report_directory_path = 'Z:\Helium_Tan\R02_PTMDIA\Pro\Spectronaut\directDIA' + '\\'   #Point to directory path
platform = 'timsTOF_Pro_stringent'                                                 #Instrument platform, can assign any identifying name
workflow = 'Spectronaut_DDALibrary'                                                #Data analysis workflow used, can assign any identifying name
ptm = '[+80]'                                                                #PTM annotation, please include parentheses
#Assign PTM localization threshold. Regular = 0.75, Stringent = 0.99 for Spectronaut (Lou et al., Nature Communications, 2023)
localization_thresholds_regular = 0.01
localization_thresholds_stringent = 0.51
conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]         #Create a list of the spike level concentrations, as floats (not strings!)
spiked_peptides = lights+heavies                                                   #List of spiked peptides that need to be searched for in data, PEPTIDES MUST HAVE SAME ANNOTATIONS AS USED IN SPECTRONAUT
#End of user input section

spectronaut = pd.read_csv(report_directory_path + '20230227_104026_PhosphoDIA_R02_Report_directDIA.tsv', sep = '\t', low_memory = False)
singly_phosphorylated = [x for x in list(spiked_peptides) if str(x).count(ptm)==1]               #Keep sequences that only have one instance of the PTM
print('read')


def localize(spectronaut,singly_phosphorylated,localization_thresholds_regular,localization_thresholds_stringent):
    '''

    Args:
        spectronaut: Spectronaut dataframe
        singly_phosphorylated
        localization_thresholds_stringent
        localization_thresholds_regular

    Returns:
        Bar plot indicating percent phosphosites retained across the titration curve for regular and stringent confidence thresholds.
        Kernel density plot showing distribution of PTM confidence scores from 0-1 across the titration curve

    '''

    dist = {}
    localized = {}

    for c in conditions:

        probabilities = []

        conc = spectronaut.loc[spectronaut['R.Condition'].str.contains(str(c))]  #Only runs from that spiked condition
        spiked = conc.loc[conc['FG.IntMID'].isin(singly_phosphorylated)]         #Only spiked peptide entries

        for index, row in spiked.iterrows():

            probs = row['EG.PTMProbabilities [Phospho (STY)]']                   #Phospho probabilities, will need to select appropriate row for different PTM

            if probs != float:
                prob = probs.split(';')

            seq = row['PEP.StrippedSequence']       #Stripped sequence
            mod_seq = row['EG.IntPIMID']            #Modified sequence

            # Finding probabilitiy info

            all_sty = []
            find_sty = re.finditer('[STY](\[\+80\])?',mod_seq)  #Find position of S|T|Y, regardless of if followed by +80 or not, this regex needs to change depending on PTM
            for match in find_sty:
                all_sty.append(match.start())

            all_phos = []
            find_phos = re.finditer('[STY](\[\+80\])',mod_seq)  #Find position of phosphorylated S|T|Y, which are followed by +80, this regex needs to change depending on PTM
            for match in find_phos:
                all_phos.append(match.start())

            index = [all_sty.index(x) for x in all_phos][0]
            actual_prob = prob[index]                           #Probability for actually modified amino acid residue

            probabilities.append(float(actual_prob))

        reg = [x for x in probabilities if x >= localization_thresholds_regular]                    #All probabilities that meet regular PTM site confidence threshold
        stringent = [x for x in probabilities if x >= localization_thresholds_stringent]            #All probabilities that meet stringent PTM site confidence threshold

        perc_reg = (len(reg) / len(probabilities)) * 100                                            #Percent entries over regular threshold
        perc_stringent = (len(stringent) / len(probabilities)) * 100                                #Percent entries over stringent threshold
        print(len(probabilities))

        localized[c] = (perc_reg, perc_stringent)

        dist[c] = probabilities

    df = pd.DataFrame.from_dict({key: pd.Series(value) for key, value in dist.items()})
    localized_df = pd.DataFrame(localized, index=['Regular', 'Stringent'])
    localized_df = localized_df.transpose()
    localized_df.to_csv(report_directory_path + 'SiteLocalization.tsv', sep='\t', index=False)

    ax = df.plot.kde()                       #Bar chart indicating percent of phosphosites retained at different thresholds
    bx = localized_df.plot.bar()             #Kernel density plots for confidence score distributions across the titration curve

    plt.show()


localize(spectronaut,singly_phosphorylated,localization_thresholds_regular,localization_thresholds_stringent)
