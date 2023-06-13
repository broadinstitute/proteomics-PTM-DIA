import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os
import matplotlib.pyplot as plt



#For TV, but can be used to read in a file with spiked peptides as a column
lights = list(pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Lights.tsv', delimiter= '\t')['DIANN_Mods'])        #Modified sequence annotates phosphopeptides as they appear on DIA-NN
heavies = list(pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Heavies.tsv', delimiter= '\t')['DIANN_Mods'][0:233])

#USER INPUT: Assign variables here#
report_directory_path = 'Z:\Helium_Tan\R02_PTMDIA\Pro\DIANN\Predicted_Library\dia_nn\out' + '\\'    #Point to directory path
platform = 'timsTOF_Pro'                                                 #Instrument platform, can assign any identifying name
workflow = 'DIANN_DDALibrary'                                                #Data analysis workflow used, can assign any identifying name
ptm = '(UniMod:21)'                                                                #PTM annotation, please include parentheses
#Assign PTM localization threshold. Regular = 0.01, Stringent = 0.51 for DIA-NN (Lou et al., Nature Communications, 2023)
localization_thresholds_regular = 0.75
localization_thresholds_stringent = 0.99
conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]         #Create a list of the spike level concentrations, as floats (not strings!)
spiked_peptides = lights+heavies                                                   #List of spiked peptides that need to be searched for in data, PEPTIDES MUST HAVE SAME ANNOTATIONS AS USED IN SPECTRONAUT
#End of user input section

diann = pd.read_csv(report_directory_path + 'report.tsv', delimiter='\t', low_memory=False)      #Read in report
singly_phosphorylated = [x for x in list(spiked_peptides) if str(x).count(ptm)==1]               #Keep sequences that only have one instance of the PTM
print('read')


#USER WILL LIKELY HAVE TO ADAPT THIS SECTION
def convert_diann(diann_df):
    '''

    Args:
        diann_df: DIA-NN report dataframe

    Returns: DIA-NN data frame with a new column that indicates spike level
    '''
    spike = []

    for r in diann_df['Run']:
        parts = r.split('_')

        fmol = [i for i in parts if 'fmol' in i]

        if len(fmol) == 0:  # NoSpikeRuns
            amt = 0

        else:
            fetch = fmol[0].split('fmol')[0]
            amt = fetch.replace('p', '.')

        spike.append(float(amt))

    diann_df['Spike'] = spike                   #Add additional column to the diann report
    return(diann_df)

#Update DIA-NN dataframe
diann = convert_diann(diann)


def localize(diann,singly_phosphorylated,localization_thresholds_regular,localization_thresholds_stringent):
    '''

    Args:
        diann: DIA-NN dataframe with new 'Spike' column
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
        conc = diann.loc[diann['Spike'] == c]                                       #Only runs from that spiked condition
        spiked = conc.loc[conc['Modified.Sequence'].isin(singly_phosphorylated)]    #Only spiked peptide entries
        probabilities = list(spiked['PTM.Site.Confidence'])

        dist[c] = probabilities

        reg = [x for x in probabilities if x >= localization_thresholds_regular]            #All probabilities that meet regular PTM site confidence threshold
        stringent = [x for x in probabilities if x >= localization_thresholds_stringent]    #All probabilities that meet stringent PTM site confidence threshold

        perc_reg = (len(reg) / len(probabilities)) * 100                                    #Percent entries over regular threshold
        perc_stringent = (len(stringent) / len(probabilities)) * 100                        #Percent entries over stringent threshold

        localized[c] = (perc_reg, perc_stringent)
        print(c, len(probabilities), perc_reg, perc_stringent)

    df = pd.DataFrame.from_dict({key: pd.Series(value) for key, value in dist.items()})
    localized_df = pd.DataFrame(localized, index=['Regular', 'Stringent'])
    localized_df = localized_df.transpose()
    localized_df.to_csv(report_directory_path + 'SiteLocalization.tsv', sep='\t')

    ax = df.plot.kde()                      #Bar chart indicating percent of phosphosites retained at different thresholds
    bx = localized_df.plot.bar()            #Kernel density plots

    plt.show()

localize()
