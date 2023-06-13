import pandas as pd
import matplotlib
from matplotlib_venn import venn2, venn3
import matplotlib.pyplot as plt
from upsetplot import plot
import seaborn as sns
import numpy as np
import math

#For TV, but can be used to read in a file with spiked peptides as a column
lights = list(pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Lights.tsv', delimiter= '\t')['DIANN_Mods'])        #Modified sequence annotates phosphopeptides as they appear on DIA-NN
heavies = list(pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Heavies.tsv', delimiter= '\t')['DIANN_Mods'][0:233])

#USER INPUT: Assign variables here#
report_directory_path = 'Z:\Helium_Tan\R02_PTMDIA\Pro\DIANN\DDA_Library\Three_varmods\dia_nn\out' + '\\'    #Point to directory path
platform = 'timsTOF_Pro'                                                                                    #Instrument platform, can assign any identifying name
workflow = 'DIANN_DDALibrary'                                                                               #Data analysis workflow used, can assign any identifying name
spiked_peptides = lights+heavies
conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]         #Create a list of the spike level concentrations, as floats (not strings!)
diann = pd.read_csv(report_directory_path + 'report.tsv', sep = '\t')


#USER WILL LIKELY HAVE TO ADAPT THIS SECTION
def convert_diann(diann_df):
    '''

    Args:
        diann_df: DIA-NN report dataframe

    Returns: DIA-NN data frame with a new column that indicates spike level
    '''
    spike = []

    for r in diann['Run']:
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
spiked = diann.loc[diann['Modified.Sequence'].isin(spiked_peptides)]            #Keep rows from spiked peptides


def characterize_lost(spiked, conditions):
    '''
    Args:
        spiked: subset of DIA-NN that is for spiked peptides
        conditions: list of spike levels

    Returns: Kernel density plot of any parameter for peptides that are maintained vs lost across the titration curve
    '''

    highest = spiked.loc[spiked['Spike'] == conditions[-1]]       #All entries from highest spike level
    lowest = spiked.loc[spiked['Spike'] == conditions[0]]         #All entries from lowest spike level

    #What's lost along the way?
    unique = set(highest['Modified.Sequence']) - set(lowest['Modified.Sequence'])       #Sequences that are at the highest concentration and not in the lowest (unique to highest concentration)

    highest_only = highest.loc[highest['Modified.Sequence'].isin(unique)]               #Subset of report with sequences that are unique to highest concentration
    highest_shared = highest.loc[~highest['Modified.Sequence'].isin(unique)]            #Subset of report with sequences that are shared between highest and lowest concentration

    list = [highest_shared, highest_only]
    names = ['Kept', 'Lost']

    dict = {}

    for i in range(0, len(list)):
        print(i)
        df = list[i]
        name = names[i]

        df = df.dropna(subset=['Precursor.Quantity'])

        sequence = df['Modified.Sequence']
        charge = df['Precursor.Charge']
        im = df['IM']
        localization = df['PTM.Site.Confidence']
        score = df['CScore']
        quant = df['Precursor.Quantity']
        log_quant = np.log10(quant)
        cleaned_log = [v for v in log_quant if not math.isnan(v) and not math.isinf(v)]
        rt = df['RT']

        #USER MUST INDICATE WHICH OF THE ABOVE PARAMETERS NEED TO BE ASSESSED
        dict[name] = rt

    df = pd.DataFrame.from_dict({key: pd.Series(value) for key, value in dict.items()})

    ax = df.plot.kde()
    plt.show()                  #Save figure to stop running script


characterize_lost(spiked, conditions)
