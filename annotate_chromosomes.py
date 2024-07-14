#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Annotate snps (no flipping)
"""

import pandas as pd 
import numpy as np 
import argparse
import sys
import warnings
warnings.filterwarnings('ignore')

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='adds RSIDs etc from a key file' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='input file',required='True')
    parser.add_argument('-k', '--key', help='RSID folder with files named chr*.ukbb.rsid.key, make sure it ends with / ', required='True')
    parser.add_argument('-o', '--output', help='output file', required='True')
    results = parser.parse_args(args)
    return (results.input , results.key, results.output )

    
def clean_me_plz(input_file, key_file , output):
    raw = pd.read_csv(input_file , delimiter=r"\s+" )
    chr_name = raw.CHROM.mode().values[0]
    key_name = key_file + "chr" + str(chr_name) + ".ukbb.rsid.key"
    key = pd.read_csv(key_name, sep="\t")
    merged = raw.merge(key,on="ID")
    merged = merged[["CHROM", "GENPOS","ID","RSID", "ALLELE1","ALLELE0","BETA","SE"]].copy()
    merged.columns = ["CHR", "POS","ID","RSID", "A1","A0","BETA","SE"]
    
    merged["Z"] = merged.BETA / merged.SE
    merged["P"] = 1 / abs(merged.Z)   
    
    merged.to_csv(output, sep="\t" , index=None )
    
if __name__ == '__main__':

    input_file, key_file , output_file = check_arg(sys.argv[1:])
    clean_me_plz(input_file, key_file , output_file)