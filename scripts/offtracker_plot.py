#!/usr/bin/env python
# -*- coding: utf-8 -*-

import offtracker.X_offplot as xoffplot
import pandas as pd
import argparse
import os

def main():
    parser = argparse.ArgumentParser()
    parser.description='Draw the plot of the off-targets with genomic sequences.\nIf .pdf file is too large, try to use .png file instead.'
    parser.add_argument('--result' , type=str, required=True,  help='The file of Offtracker_result_{outname}.csv' )
    parser.add_argument('--sgrna'  , type=str, required=True,  help='Not including PAM' )
    parser.add_argument('--pam'    , type=str, default='NGG',  help='PAM sequence. Default is "NGG".' )
    parser.add_argument('--output' , type=str, default='same', help='The output file. Default is Offtracker_result_{outname}.pdf')

    args = parser.parse_args()
    if args.output == 'same':
        dir_savefig = args.result.replace('.csv', '.pdf')
    else:
        dir_savefig = args.output

    outname = os.path.basename(args.result).replace('Offtracker_result_', '').replace('.csv', '')
    gRNA = args.sgrna
    PAM = args.pam
    full_seq = gRNA + PAM
    
    df_result = pd.read_csv(args.result)
    n_pos = len(df_result)

    xoffplot.offtable(df_result, full_seq, length_pam = len(PAM), col_seq='target', threshold=2,
                      title=f'{outname} ({n_pos} sites)',
                      savefig=dir_savefig)
    
    return f'The plot is saved as {dir_savefig}'

if __name__ == '__main__' :
    result = main()
    print(result)
