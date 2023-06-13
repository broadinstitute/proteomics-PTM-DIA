import pandas as pd
import matplotlib
from matplotlib_venn import venn2, venn3
import matplotlib.pyplot as plt
from upsetplot import plot
import seaborn as sns
import numpy as np
import math
import re

# Read in files

# Spectronaut
SN_DDALib = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\Spectronaut\DDALibrary\\20230308_091202_R02_PhosphoDIA_Pro_DDALibrary_Report.tsv',sep='\t', low_memory=False)
SN_directDIA = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\Spectronaut\directDIA\\20230227_104026_PhosphoDIA_R02_Report_directDIA.tsv', sep='\t',low_memory=False)
SN_HybridLib = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\Spectronaut\Hybrid\\20230313_103651_R02_PhosphoDIA_Pro_HybridLibrary_Report.tsv',sep='\t', low_memory=False)

# DIA-NN
DIANN_DDALib = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\DIANN\DDA_Library\Three_varmods\dia_nn\out\\report.tsv',sep='\t')
DIANN_CombinedLib = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\DIANN\Combined_DDA_Predicted_Library\dia_nn\out\\report.tsv', sep='\t')
DIANN_PredictedLib = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\DIANN\Predicted_Library\dia_nn\out\\report.tsv',sep='\t')


def convert(sequence):
    """
        Args: Sequence in Spectronaut annotation format.
        Returns: Sequence converted to reflect DIA-NN annotation format.

    """
    dict = {

        '_\[\+42\]': '(UniMod:1)',
        '\[\+57\]': '(UniMod:4)',
        '\[\+80\]': '(UniMod:21)',
        '\[\+16\]': '(UniMod:35)',
        '\[\+8\]': '(UniMod:259)',
        '\[\+10\]': '(UniMod:267)'
    }

    for x in dict.keys():
        sequence = re.sub(x, dict[x], sequence)  # Update sequence with replaced annotation

    return (sequence)


# Spectronaut
SN = [SN_DDALib, SN_directDIA, SN_HybridLib]  # The Spectronaut data frames
for search in SN:
    search['Modified.Sequence'] = search['FG.IntMID'].apply(lambda x: convert(
        x))  # Apply convert function to peptide sequence column in, create new column that matches DIA-NN


def filter_sn(df):
    """
        Args: Spectronaut report.
        Returns: Stringently ocalized phosphopeptide-only report.

    """
    phospho = df.loc[df['Modified.Sequence'].str.contains('\(UniMod:21\)')]  # Keep only phosphopeptides
    filtered = phospho.loc[phospho[
                               'EG.PTMAssayProbability'] >= 0.99]  # Keep only stringently localized phosphosites (0.99 selected as value for stringent localization)
    return (filtered)


# Filtered phospho-only Spectronaut datasets
SN_DDALib = filter_sn(SN_DDALib)
SN_directDIA = filter_sn(SN_directDIA)
SN_HybridLib = filter_sn(SN_HybridLib)


# DIA-NN
def filter_diann(df):
    """
        Args: DIA-NN report.
        Returns: Stringently ocalized phosphopeptide-only report.

    """
    phospho = df.loc[df['Modified.Sequence'].str.contains('\(UniMod:21\)')]  # Keep only phosphopeptides
    filtered = phospho.loc[phospho[
                               'PTM.Site.Confidence'] >= 0.51]  # Keep only stringently localized phosphosites (0.51 selected as value for stringent localization)
    return (filtered)


# Filtered phospho-only DIA-NN datasets
DIANN_DDALib = filter_diann(DIANN_DDALib)
DIANN_CombinedLib = filter_diann(DIANN_CombinedLib)
DIANN_PredictedLib = filter_diann(DIANN_PredictedLib)

all = [DIANN_DDALib, DIANN_PredictedLib, DIANN_CombinedLib, SN_DDALib, SN_directDIA,SN_HybridLib]  # List of data frames that should be included in dataframe
names = ['DIANN_DDALib', 'DIANN_PredictedLib', 'DIANN_HybridLib', 'SN_DDALib', 'SN_directDIA','SN_HybridLib']  # Names you want for each data set in the venn diagram

sequences_all = [set(x['Modified.Sequence']) for x in all]  # Select the Modified.Sequence column from each dataframe
venn3(sequences_all, names, set_colors=('#0000ff', '#ff4500', '#006400'))  # Crate venn diagram for up to 3 dataframes


def find_unique(all, sequences_all):
    """
        Returns list where each item is a list of sequences unique to the respective dataframe.
        Args:
            all: List of dataframes
            sequences_all: List of sequences in those data frames
        Returns: List of dataframes that contain only sequences unique to each individual data frame

    """
    updated = []

    for i in range(0, len(sequences_all)):  # For each dataframe
        seqs = sequences_all[i]

        new = list(sequences_all)  # Make copy of full list
        new.remove(seqs)  # Remove the sequences that are in this dataframe

        other = set.union(*new)  # The unique sequences in the other dataframes
        unique = seqs - other  # Sequences unique to this particular dataframe
        print(len(unique))

        # keep only sequences that are unique in the current dataframe

        df = all[i]
        find_unique = df.loc[df['Modified.Sequence'].isin(unique)]  # Want full dataframe so that you have all the information relating to these unique sequences

        updated.append(find_unique)  # Add dataframe of unique sequences to the list initialized above

    return (updated)


unique = find_unique(all, sequences_all)


def compare(all, names):
    '''
    Args:
        all: List of dataframes
        names: List of names of above dataframes

    Returns: Density plot distribution of assigned parameter across different dataframes
    '''
    dict = {}
    for i in range(0, len(all)):  # Iterate through dataframes, extract information depending on the software it is from
        name = names[i]
        df = all[i]

        if 'SN' in name:  # Spectronaut
            df = df.dropna(subset=['EG.PTMAssayProbability'])

            charge = df['FG.Charge']
            im = df['FG.IonMobility']
            localization = df['EG.PTMAssayProbability']
            score = df['EG.Cscore']
            quant = df['FG.Quantity']
            log_quant = np.log10(quant)
            cleaned_log = [v for v in log_quant if not math.isnan(v) and not math.isinf(v)]
            rt = df['EG.ApexRT']

        if 'DIANN' in name:  # DIA-NN
            df = df.dropna(subset=['Precursor.Quantity'])

            charge = df['Precursor.Charge']
            im = df['IM']
            localization = df['PTM.Site.Confidence']
            score = df['CScore']
            quant = df['Precursor.Quantity']
            log_quant = np.log10(quant)
            cleaned_log = [v for v in log_quant if not math.isnan(v) and not math.isinf(v)]
            rt = df['RT']

        dict[name] = cleaned_log  # Here, need to determine which variable you would need plotted

    df = pd.DataFrame.from_dict({key: pd.Series(value) for key, value in dict.items()})
    ax = df.plot.kde()

    plt.show()

#Compare values unique to each dataframe to observe if there is any distingushing pattern
compare(unique, names)