import pandas as pd
import matplotlib
from matplotlib_venn import venn2, venn3
import matplotlib.pyplot as plt
from upsetplot import plot
import seaborn as sns
import numpy as np
import math

#For TV, but can be used to read in a file with spiked peptides as a column
lights = list(pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Lights.tsv', delimiter= '\t')['Modified'])        #Modified sequence annotates phosphopeptides as they appear on Spectronaut
heavies = list(pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Heavies.tsv', delimiter= '\t')['LibraryHeavies'][0:233])

#USER INPUT: Assign variables here#
report_directory_path = 'Z:\Helium_Tan\R02_PTMDIA\Pro\Spectronaut\DDALibrary' + '\\'    #Point to directory path
file_name = '20230308_091202_R02_PhosphoDIA_Pro_DDALibrary_Report.tsv'
platform = 'timsTOF_Pro'                                                                                    #Instrument platform, can assign any identifying name
workflow = 'DIANN_DDALibrary'                                                                               #Data analysis workflow used, can assign any identifying name
spiked_peptides = lights+heavies
conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]         #Create a list of the spike level concentrations, as floats (not strings!)
spectronaut = pd.read_csv(report_directory_path + file_name, sep = '\t')



spiked = spectronaut.loc[spectronaut['FG.IntMID'].isin(spiked_peptides)]            #Keep rows from spiked peptides
print(len(spiked))

def characterize_lost(spiked, conditions):
    '''
    Args:
        spiked: subset of Spectronaut report that is for spiked peptides
        conditions: list of spike levels

    Returns: Kernel density plot of any parameter for peptides that are maintained vs lost across the titration curve
    '''

    highest = spiked.loc[spiked['R.Condition'].str.contains(str(conditions[-1]))]       #All entries from highest spike level
    lowest = spiked.loc[spiked['R.Condition'].str.contains(str(conditions[0]))]         #All entries from lowest spike level
    print(len(highest), len(lowest))


    #What's lost along the way?
    unique = set(highest['FG.IntMID']) - set(lowest['FG.IntMID'])       #Sequences that are at the highest concentration and not in the lowest (unique to highest concentration)

    highest_only = highest.loc[highest['FG.IntMID'].isin(unique)]               #Subset of report with sequences that are unique to highest concentration
    highest_shared = highest.loc[~highest['FG.IntMID'].isin(unique)]            #Subset of report with sequences that are shared between highest and lowest concentration

    list = [highest_shared, highest_only]
    names = ['Kept', 'Lost']

    dict = {}

    for i in range(0, len(list)):
        print(i)
        df = list[i]
        name = names[i]

        df = df.dropna(subset=['EG.PTMAssayProbability'])

        charge = df['FG.Charge']
        im = df['FG.IonMobility']
        localization = df['EG.PTMAssayProbability']
        score = df['EG.Cscore']
        quant = df['FG.Quantity']
        log_quant = np.log10(quant)
        cleaned_log = [v for v in log_quant if not math.isnan(v) and not math.isinf(v)]
        rt = df['EG.ApexRT']

        #USER MUST INDICATE WHICH OF THE ABOVE PARAMETERS NEED TO BE ASSESSED
        dict[name] = rt

    df = pd.DataFrame.from_dict({key: pd.Series(value) for key, value in dict.items()})

    ax = df.plot.kde()
    plt.show()                  #Save figure to stop running script


characterize_lost(spiked, conditions)
