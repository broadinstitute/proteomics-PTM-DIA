import pandas
import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os

lights = list(pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Lights.tsv', delimiter= '\t')['Modified'])       #Modified sequence annotates phosphopeptides as they appear on Spectronaut
heavies = list(pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Heavies.tsv', delimiter= '\t')['LibraryHeavies'][0:233])

#USER INPUT: Assign variables here

report_directory_path = 'Z:\Helium_Tan\R02_PTMDIA\Pro\Spectronaut\DDALibrary\\'    #Point to directory path
file_name = '20230308_091202_R02_PhosphoDIA_Pro_DDALibrary_Report.tsv'             #File name in tsv format, R.Condition column must be spike level concentration, e.g. '10.0' for 10 fmol
platform = 'timsTOF_Pro_stringent'                                                 #Instrument platform, can assign any identifying name
workflow = 'Spectronaut_DDALibrary'                                                #Data analysis workflow used, can assign any identifying name
localization_threshold = 0.99                                                      #Assign PTM localization threshold. Regular = 0.75, Stringent = 0.99 (Lou et al., Nature Communications, 2023)
conditions = [0.004, 0.01, 0.02, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]         #Create a list of the spike level concentrations, as floats (not strings!)
norm_spike = 1.0                                                                   #The spike level concentration that quantities at other spike level quantities will be ratio'd to
spiked_peptides = lights+heavies                                                   #List of spiked peptides that need to be searched for in data, PEPTIDES MUST HAVE SAME ANNOTATIONS AS USED IN SPECTRONAUT

#Read in report and filter for PTM site localization
spectronaut = pd.read_csv(report_directory_path + file_name, delimiter='\t', low_memory=False)      #Read in Spectronaut report
spectronaut = spectronaut[spectronaut['EG.PTMAssayProbability'] >= localization_threshold]          #Filter data set for entries that meet the PTM site localization confidence score threshold
print('read')

def makePath(filepath):
    '''
    Makes filepath
    '''
    isExist = os.path.exists(filepath)
    if not isExist:
        os.makedirs(filepath)

def summarize_quant(spectronaut, spiked_peptides, conditions):
    summary = {}

    for sequence in spiked_peptides:

        find = sequence  # Spiked in peptide sequence
        found_at_spike = 0

        single = spectronaut.loc[spectronaut['FG.IntMID'] == find]  # Only look at the entries where that spiked peptide was found
        single = single[['R.Condition', 'R.Replicate', 'FG.IntMID','FG.Quantity']]  # Keep only the column indicating condition/spike, peptide sequence, and quantity
        single = single[single['FG.Quantity'].notnull() & single['FG.Quantity'] > 0]  # Remove columns where quantity is null and quantity is greater than 0

        # Add new column with spiked amount: Need to check this one out
        spiked = []
        for x in single['R.Condition']:
            amount = float(x.replace('fmol', ''))
            spiked.append(amount)

        single['Spiked'] = spiked

        df_list = []
        for x in conditions:  # For every spike level on the curve, create a new data frame, iterate through each spike level

            spike_label = str(x)
            if spike_label not in summary:
                summary[spike_label] = {}
                summary[spike_label]['Reps'] = []
                summary[spike_label]['CVs'] = []

            data = {}
            current = single[single['Spiked'] == x]  # New df where you are only looking at entries associated with that spiked amount

            # Gather information
            data['Instrument'] = platform
            data['Workflow'] = workflow
            data['Peptide'] = find
            data['Spike'] = x
            data['Log Spike'] = math.log10(x)  # Plot log spiked amount vs log quantity

            # If spiked peptide is not found at this point on the curve
            if len(current) == 0:
                data['Quant_EachRep'] = []
                data['Number of Reps'] = 0
                data['Mean Quant'] = None
                data['Log Mean Quant'] = None
                data['Percent CV'] = None
                data['Expected Ratio'] = norm_spike / x

            # If spiked peptide is found at this point on the curve
            if len(current) > 0:
                found_at_spike += 1
                grouped = current.groupby('R.Replicate', as_index=False).agg({'R.Condition': 'first', 'Spiked': 'first', 'FG.IntMID': 'first','FG.Quantity': 'sum'})  # Condense to one line per replicate, sum quantities
                num_reps = len(grouped['R.Replicate'])  # Number of replicates it's found in
                quant = math.log10(mean(grouped['FG.Quantity']))  # Take the log of the mean quantity
                cv = (grouped['FG.Quantity'].std() / mean(grouped['FG.Quantity'])) * 100  # Calculate percent CV between quantities found across replicates

                # Add relevant data to dictionary
                data['Quant_EachRep'] = [x for x in grouped['FG.Quantity']]
                data['Number of Reps'] = num_reps
                data['Mean Quant'] = mean(grouped['FG.Quantity'])
                data['Log Mean Quant'] = quant
                data['Percent CV'] = cv
                data['Expected Ratio'] = norm_spike / x

                # Add reps and CVs information to summarizing table
                summary[spike_label]['Reps'].append(num_reps)
                summary[spike_label]['CVs'].append(cv)

                df_list.append(data)  # Collect all dictionaries into a list, so by the end it will be a list of dicts for each point on the curve

        df = pd.DataFrame(df_list)  # Combine all dicts into one data frame and export separately for each spiked peptide
        if found_at_spike > 0:
            info_path = report_directory_path + platform + '_AllSpikes_Info/'
            makePath(info_path)
            x = find.replace(':', '')
            df.to_csv(info_path + x + workflow + '_output.tsv', sep='\t')

        # Ratios: This is done where each injection a peptide is found in, the quant is ratio'd to the mean quant at norm_spike
        if len(df) != 0:
            ratios = []
            set = df[df['Spike'] == norm_spike]  # Find information about the spiked peptide at the set ratio level (norm_spike), using data frame created above

            if len(set) > 0:  # If there is quant at norm_spike level
                set_mean = mean(
                    set['Quant_EachRep'].tolist()[0])  # Get the mean quant of peptide at the norm_spike level

                for index, row in df.iterrows():  # Add actual ratios of norm_spike to peptide as new column to dataframe
                    current = row['Quant_EachRep']
                    if isinstance(current, float) == False:
                        rat = [set_mean / x for x in current]
                        ratios.append(rat)
                    else:
                        ratios.append(None)

                df['Actual Ratios'] = ratios

                # Create new dataframe that will be used to plot replicate ratios to norm_spike for each spiked peptide, each row is the peptide found in any replicate at any spike level
                new_df = pd.DataFrame(columns=['Peptide', 'Spike', 'Expected Ratio', 'Actual Ratio'])
                for index, row in df.iterrows():
                    peptide = row['Peptide']
                    spike = row['Spike']
                    expected = row['Expected Ratio']

                    if row['Actual Ratios'] != None:
                        for x in row['Actual Ratios']:
                            sub = {'Peptide': peptide, 'Spike': spike, 'Expected Ratio': expected, 'Actual Ratio': x}
                            sub = pd.DataFrame.from_dict([sub])
                            new_df = pandas.concat([new_df, sub])

                all_spikes_path = report_directory_path + platform + '_AllSpikes_Replicates/'
                makePath(all_spikes_path)
                new_df.to_csv(all_spikes_path + find + '_replicates_output.tsv',sep='\t')  # Export second set of dataframes that can be used for creating the actual boxplots

summarize_quant(spectronaut, spiked_peptides, conditions)