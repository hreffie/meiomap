#!/usr/bin/env python3
import numpy as np
import pandas as pd
import argparse
import os


# Arguement parser
def parse_arguments():
    parser = argparse.ArgumentParser(
        description='''Simple implementation of the pilot stage of haplotype phasing of oocytes. The input file is a tab
           delimited text file.''')
    parser.add_argument('--infile', dest='infile', type=argparse.FileType('r', encoding='UTF-8'),
                        required=True,
                        help="Input file is a tab delimited text file containing Chr, Pos, Maternal Genotype, Reference Genotype, Cell Genotype")
    parser.add_argument('--outdir', dest='outdir', required=True,
                        help="Output directory. Default: cwd")
    args = parser.parse_args()
    infile = args.infile
    outdir = args.outdir
    return infile, outdir


# Preprocessing of data
def preproc(infile):
    # View all columns
    pd.set_option('display.max_columns', None)

    # Read input data
    df = pd.read_csv(infile, sep="\t", names=['Name', 'Chr', 'Position', 'gDNAGType', '1PB1GType', '1PB2GType',
                                              '1eggGType', '2PB1GType', '2PB2GType', '2eggGType',
                                              '3PB1GType', '3PB2GType', '3eggGType', '4PB1GType', '4PB2GType',
                                              '4eggGType'])

    # Extract relevant maternal genotype (heterozygous mothers)
    data = df[df['gDNAGType'] == 'AB']

    # Extract reference genotype
    ref = data['1eggGType']

    # Reorder reference columns for clarity
    data.insert(4, 'Ref', ref)

    # Change GT to numbers and filter out NCs
    data.replace(to_replace='AA', value='0', inplace=True)
    data.replace(to_replace='AB', value='1', inplace=True)
    data.replace(to_replace='BB', value='2', inplace=True)
    data.replace('NC', np.NaN, inplace=True)
    data.dropna(inplace=True)

    # Create a list of empty list to be filled in by chr,pos and phase.

    cellnames = '_1PB1', '_1PB2', 'Egg1', '_2PB1', '_2PB2', 'Egg2', '_3PB1', '_3PB2', 'Egg3', '_4PB1', '_4PB2', 'Egg4'
    _1PB1 = []
    _1PB2 = []
    Egg1 = []
    _2PB1 = []
    _2PB2 = []
    Egg2 = []
    _3PB1 = []
    _3PB2 = []
    Egg3 = []
    _4PB1 = []
    _4PB2 = []
    Egg4 = []
    cells = [_1PB1, _1PB2, Egg1, _2PB1, _2PB2, Egg2, _3PB1, _3PB2, Egg3, _4PB1, _4PB2, Egg4]
    x = 0
    # Compare SNPs to REF for phasing and add to empty list of each cell
    for x in range(0, 11):
        for index, row in data.iterrows():
            if row.iloc[5 + x] == row.iloc[4]:
                cells[x].append('1')
            elif row.iloc[5 + x] == '1' and (row.iloc[4] == '2' or row.iloc[4] == '0'):
                cells[x].append('0.5')
            elif (row.iloc[5 + x] == '2' or row.iloc[5 + x] == '0') and row.iloc[4] == '1':
                cells[x].append('0.5')
            else:
                cells[x].append('0')

        # Merge data frames and set stop pos = start pos
        cells[x] = pd.DataFrame(
            {'Chr': data['Chr'], 'Start': data['Position'], 'Stop': data['Position'], 'Phase': cells[x]})
        cells[x].reset_index(drop=True, inplace=True)
    return cells


# Cluster range of phases
def cluster_phase(infile):
    cells = preproc(infile)
    # Cluster phase areas
    # _1PB1, _1PB2, Egg1, _2PB1, _2PB2, Egg2, _3PB1, _3PB2, Egg3, _4PB1, _4PB2, Egg4 = [i for i in cells]
    #
    # Set first start position in phase equal to all start position in phase for phase range
    for x in range(0, 11):
        for i in range(0, len(cells[x]) - 1):
            if cells[x].loc[i, 'Phase'] == cells[x].loc[i + 1, 'Phase']:
                cells[x].loc[i + 1, 'Start'] = cells[x].loc[i, 'Start']
                # Remove duplicates
                cells[x].drop_duplicates(subset=['Start'], keep='last')

    return cells


def write_bed(infile, outdir):
    files = cluster_phase(infile)
    filenames = ['_1PB1', '_1PB2', 'Egg1', '_2PB1', '_2PB2', 'Egg2', '_3PB1', '_3PB2', 'Egg3', '_4PB1', '_4PB2', 'Egg4']

    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print("Output directory will be created ", outdir)

    for x in range(0, len(files) - 1):
        temp = os.path.join(outdir, '{}.bed'.format(filenames[x]))
        files[x].to_csv(temp, sep='\t')


def main():
    infile, outdir = parse_arguments()
    write_bed(infile, outdir)


if __name__ == '__main__':
    main()
