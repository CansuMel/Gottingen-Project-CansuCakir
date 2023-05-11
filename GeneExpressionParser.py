#!/usr/bin/env python3
"""
Date last edited 2023-04-21 19:08
Description: Expression Analysis of Potential IFNβ Bioactivity Biomarkers in Patients with Multiple Sclerosis
@Author: Cansu M. Cakir
"""
# 1.Parse the command line
import argparse

parser = argparse.ArgumentParser(
    prog="GeneExpressionParser",
    description="Prints expression patterns of genes from microarray data upon treatment with IFB1a & IFB1b based on Log2FC value",
    epilog="If there are problems contact Mel!"
)
parser.add_argument("-f", '--filename', required=True,
                    help='Only reads .tsv / tab-separated file with column name "Factor Value[compound]"!')
args = parser.parse_args()

# 2.Open the file and read the contents
with open(args.filename, "r") as file:
    content = file.readlines()

# 3.Parse all the content into a dataframe for experiment design
import pandas as pd

ExpDesign_df = pd.read_csv(args.filename, sep='\t')

# 4.Group Treatment 1, Treatment 2 and Control
IFB1a = ExpDesign_df[ExpDesign_df["Factor Value[compound]"].str.contains("Interferon beta-1a") == True]
IFB1b = ExpDesign_df[ExpDesign_df["Factor Value[compound]"].str.contains("Interferon beta-1b") == True]
control = ExpDesign_df[ExpDesign_df["Factor Value[compound]"].str.contains("none") == True]

# 5.Match the control with corresponding treatment in a df based on patient number
# Write the full file name in df
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
import numpy as np


def log2fc_calculator(df):
    result_df = pd.DataFrame()
    for idx, row in df.iloc[0:].iterrows():
        # Call for the treatment and control data files
        treatment_file = df["Treatment"].iloc[idx]
        control_file = df["Control"].iloc[idx]
        treatment_data = pd.read_csv(treatment_file, sep='\t', header=None, names=['Gene_ID', 'Value'])
        treatment_data = treatment_data[1:]
        control_data = pd.read_csv(control_file, sep='\t', header=None, names=['Gene_ID', 'Value'])
        control_data = control_data[1:]
        # Merge the treatment and control data on Gene Ref_ID
        merged_data = pd.merge(treatment_data, control_data, on='Gene_ID')
        # Convert 'Value_x' and 'Value_y' columns to numeric
        merged_data['Value_x'] = pd.to_numeric(merged_data['Value_x'])
        merged_data['Value_y'] = pd.to_numeric(merged_data['Value_y'])
        # Calculate the fold change
        merged_data['Fold_Change'] = merged_data['Value_x'] / merged_data['Value_y']
        # Add the fold change column to the result dataframe
        result_df['Gene_ID'] = merged_data['Gene_ID']
        result_df[f'Fold_Change{idx}'] = merged_data['Fold_Change']

    # Calculate and append to a new column in result_df the mean of foldchange results and the log2FC of mean
    result_df['Mean_FC'] = result_df.mean(axis=1)
    result_df['Log2_FC'] = np.log2(result_df['Mean_FC'])
    return result_df


# Name the log2FC results of treatments
result_df_1a = log2fc_calculator(comp_df_1a)
result_df_1b = log2fc_calculator(comp_df_1b)

# 7.Write a result on screen showing the expression pattern and how many genes show this trend.

# Find online and replace Ref num with actual gene name XXXunsuccesfulXXX
# gene_name.replace("AFFX - BioB - 5_at", "")..

# Create empty lists for necessary expression patterns
up_regulated_both = []
up_regulated_1a = []
up_regulated_1b = []
down_regulated_both = []
none = []

# Loop through each row in result_df_1a and result_df_1b simultaneously
for idx, row_1a in result_df_1a.iterrows():
    row_1b = result_df_1b.iloc[idx]

    # Determine and define the expression pattern for iterating rows in each dataframe
    if row_1a['Log2_FC'] > 1:
        pattern_1a = 'up'
    elif row_1a['Log2_FC'] < -1:
        pattern_1a = 'down'
    else:
        pattern_1a = 'none'

    if row_1b['Log2_FC'] > 1:
        pattern_1b = 'up'
    elif row_1b['Log2_FC'] < -1:
        pattern_1b = 'down'
    else:
        pattern_1b = 'none'

    # Add the gene to the appropriate list based on the expression pattern
    if pattern_1a == 'up' and pattern_1b == 'up':
        up_regulated_both.append(row_1a['Gene_ID'])
    elif pattern_1a == 'up' and pattern_1b == 'down':
        up_regulated_1a.append(row_1a['Gene_ID'])
    elif pattern_1a == 'down' and pattern_1b == 'up':
        up_regulated_1b.append(row_1b['Gene_ID'])
    elif pattern_1a == 'down' and pattern_1b == 'down':
        down_regulated_both.append(row_1a['Gene_ID'])
    else:
        none.append(row_1a['Gene_ID'])

# Print the results
print(f"Expression pattern 1: Upregulated for IFB1a & IFB1b, {len(up_regulated_both)} genes")
print(f"Expression pattern 2: Upregulated for IFB1a, Downregulated for IFB1b, {len(up_regulated_1a)} genes")
print(f"Expression Pattern 3: Downregulated for IFB1a, Upregulated for IFB1b, {len(up_regulated_1b)} genes")
print(f"Expression Pattern 4: Downregulated for IFB1a & IFB1b, {len(down_regulated_both)} genes")

exit(0)
