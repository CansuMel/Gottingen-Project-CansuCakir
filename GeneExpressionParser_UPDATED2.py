#!/usr/bin/env python3
"""
Date last edited 2023-06-01 16:45
Description: Expression Analysis of Potential IFNβ Bioactivity Biomarkers in Patients with Multiple Sclerosis
@Author: Cansu M. Cakir
"""
# 1.Parse the command line
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(
    prog="GeneExpressionParser",
    description="Prints expression patterns of genes from microarray data upon treatment with IFB1a & IFB1b based on Log2FC value",
    epilog="Only reads .tsv / tab-separated file with column name 'Factor Value[compound]'!"
)
parser.add_argument("-f", '--filename', required=True,
                    help='Only reads .tsv / tab-separated file with column name "Factor Value[compound]"!')
args = parser.parse_args()

# 2.Open the file and read the contents
with open(args.filename, "r") as file:
    content = file.readlines()

# 3.Parse all the content into a dataframe for experiment design
ExpDesign_df = pd.read_csv(args.filename, sep='\t')

# 4.Group Treatment 1, Treatment 2 and Control
IFB1a = ExpDesign_df[ExpDesign_df["Factor Value[compound]"].str.contains("Interferon beta-1a") == True]
IFB1b = ExpDesign_df[ExpDesign_df["Factor Value[compound]"].str.contains("Interferon beta-1b") == True]
control = ExpDesign_df[ExpDesign_df["Factor Value[compound]"].str.contains("none") == True]

# 5.Match the control with corresponding treatment in a df based on patient number
# Write the full file name (with extension) in the new dataframe
IFB1a_3mo = IFB1a[IFB1a["Factor Value[time]"].str.contains("3 month") == True]
merged_df1 = pd.merge(IFB1a_3mo, control, on='Sample Characteristic[individual]')
comp_df_1a = merged_df1[["Assay_x", "Assay_y"]]
comp_df_1a.columns = ['Treatment', 'Control']
comp_df_1a = comp_df_1a + "_sample_table.txt"

IFB1b_3mo = IFB1b[IFB1b["Factor Value[time]"].str.contains("3 month") == True]
merged_df2 = pd.merge(IFB1b_3mo, control, on='Sample Characteristic[individual]')
comp_df_1b = merged_df2[["Assay_x", "Assay_y"]]
comp_df_1b.columns = ['Treatment', 'Control']
comp_df_1b = comp_df_1b + "_sample_table.txt"

# 6.Write a log2FC calculator function for two Treatment dataframes

# Foldchange = Treatment / Control
# UPDATED**  "Log2(FC)" = log2(mean(Treatment)/ log2(mean(Control)))

def log2fc_calculator(df):
    # Create an empty dataframe to store the results
    result_df = pd.DataFrame()

    for idx, row in df.iterrows():
        #Call for the treatment file and open it
        treatment_file = row['Treatment']
        treatment_data = pd.read_csv(treatment_file, sep='\t', header=None, names=['Gene_ID', 'Value'])[1:]
        treatment_data['Value'] = pd.to_numeric(treatment_data['Value'], errors='ignore')
        # Assign the values in each treatment file to a new indexed column in result_df
        result_df['T' + str(idx + 1)] = treatment_data['Value']
    # Calculate the mean result of replicas for each gene and assign it to a new column named 'T_Mean'
    result_df['T_Mean'] = result_df.mean(axis=1)
    result_df.insert(0, 'Gene_ID', treatment_data['Gene_ID'])

    for idx, row in df.iterrows():
        #Call for the control file and open it
        control_file = row['Control']
        control_data = pd.read_csv(control_file, sep='\t', header=None, names=['Gene_ID', 'Value'])[1:]
        control_data['Value'] = pd.to_numeric(control_data['Value'], errors='ignore')
        # Assign the values in each control file to a new indexed column in result_df
        result_df['C' + str(idx + 1)] = control_data['Value']

    # Calculate the mean result of replicas for each gene and assign it to a new column named 'C_Mean'
    result_df['C_Mean'] = result_df[['C1', 'C2', 'C3', 'C4']].mean(axis=1)
    # Calculate the log2 of the averages and take the fold change
    result_df['Log2_x'] = np.log2(result_df['T_Mean'])
    result_df['Log2_y'] = np.log2(result_df['C_Mean'])
    result_df['Log2_FC'] = result_df['Log2_x'] / result_df['Log2_y']

    return result_df

# Name the log2FC results of IFB1a and IFB1b treatments
result_df_1a = log2fc_calculator(comp_df_1a)
result_df_1b = log2fc_calculator(comp_df_1b)

# Create empty lists for the expression patterns
up_regulated_both = []
up_regulated_1a = []
up_regulated_1b = []
down_regulated_both = []
none = []

for idx in result_df_1a.index:
    # Access the Log2_FC value of each gene in the two treatments
    log2fc_1a = result_df_1a.loc[idx, 'Log2_FC']
    log2fc_1b = result_df_1b.loc[idx, 'Log2_FC']

    # Determine and define the expression pattern for each log2FC value
    if log2fc_1a > 1:
        pattern_1a = 'up'
    elif log2fc_1a < -1:
        pattern_1a = 'down'
    else:
        pattern_1a = 'none'

    if log2fc_1b > 1:
        pattern_1b = 'up'
    elif log2fc_1b < -1:
        pattern_1b = 'down'
    else:
        pattern_1b = 'none'

    # Add the gene to the appropriate list based on the expression pattern
    if pattern_1a == 'up' and pattern_1b == 'up':
        up_regulated_both.append(result_df_1a.loc[idx, 'Gene_ID'])
    elif pattern_1a == 'up' and pattern_1b == 'down':
        up_regulated_1a.append(result_df_1a.loc[idx, 'Gene_ID'])
    elif pattern_1a == 'down' and pattern_1b == 'up':
        up_regulated_1b.append(result_df_1b.loc[idx, 'Gene_ID'])
    elif pattern_1a == 'down' and pattern_1b == 'down':
        down_regulated_both.append(result_df_1a.loc[idx, 'Gene_ID'])
    else:
        none.append(result_df_1a.loc[idx, 'Gene_ID'])

# Print the results
print(f"Expression pattern 1: Upregulated for IFB1a & IFB1b, {len(up_regulated_both)} genes")
print(f"Expression pattern 2: Upregulated for IFB1a, Downregulated for IFB1b, {len(up_regulated_1a)} genes")
print(f"Expression Pattern 3: Downregulated for IFB1a, Upregulated for IFB1b, {len(up_regulated_1b)} genes")
print(f"Expression Pattern 4: Downregulated for IFB1a & IFB1b, {len(down_regulated_both)} genes")

exit(0)


