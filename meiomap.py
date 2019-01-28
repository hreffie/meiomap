#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os


#View all columns
pd.set_option('display.max_columns', None)

#Read input data
df = pd.read_csv('input_data_trios.txt', sep="\t", names=['Name', 'Chr', 'Position', 'gDNAGType','1PB1GType', '1PB2GType',
       '1eggGType', '2PB1GType', '2PB2GType', '2eggGType',
        '3PB1GType', '3PB2GType', '3eggGType', '4PB1GType',
        '4PB2GType', '4eggGType'])

#Extract relevant maternal genotype (heterozygous mothers)
data= df[df['gDNAGType']=='AB']

#Extract reference genotype
ref = data['1eggGType']

#Reorder columns for clarity
data.insert(4,'Ref',ref)


#Change GT to numbers and filter out NCs
data.replace(to_replace='AA', value='0', inplace=True)
data.replace(to_replace='AB', value='1', inplace=True)
data.replace(to_replace='BB', value='2', inplace=True)
print(len(data))
data.replace('NC', np.NaN, inplace=True)
data.dropna(inplace=True)

#Create a list of empty list to be filled in by chr,pos and phase.
cells = [[] for e in range(12)]
cellnames = '_1PB1', '_1PB2','Egg1', '_2PB1', '_2PB2', 'Egg2','_3PB1','_3PB2', 'Egg3', '_4PB1', '_4PB2','Egg4'

i=0
#Compare SNPs to REF
for x in cells:
    for index, row in data.iterrows():
        if row.iloc[5+i]==row.iloc[4]:
            cells[i].append('1')
        elif row.iloc[5+i]=='1' and (row.iloc[4]=='2' or row.iloc[4]=='0'):
            cells[i].append('0.5')
        elif (row.iloc[5+i]=='2' or row.iloc[5+i]=='0') and row.iloc[4] == '1':
            cells[i].append('0.5')
        else:
            cells[i].append('0')
    cells[i] = pd.DataFrame({'Chr':data['Chr'],'Start':data['Position'],'Phase':cells[i]})
    i=i+1

#Cluster phase areas
_1PB1, _1PB2, Egg1, _2PB1, _2PB2, Egg2, _3PB1, _3PB2, Egg3, _4PB1, _4PB2, Egg4 = [i for i in cells]

#Set default values for Stop positions and drop indices
_1PB1.insert(2,'Stop', _1PB1['Start'])
_1PB1.reset_index(drop=True, inplace=True)

#Set first start position in phase equal to all start position in phase for phase range
for row in range(1,len(_1PB1)-1):
    if _1PB1.loc[row,'Phase']==_1PB1.loc[row+1,'Phase']:
        _1PB1.loc[row+1,'Start']= _1PB1.loc[row,'Start']
        
#Remove duplicates
result_df = _1PB1.drop_duplicates(subset=['Start'], keep='last')
print(result_df)
