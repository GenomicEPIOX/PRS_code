import os
import sys
import pandas as pd
import numpy as np
import argparse


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Convert summary stats for PRS-cs format' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='The input sumstats file',required='True')
    parser.add_argument('-o', '--output', help='output_name',required='True')
    parser.add_argument('-a1', '--a1', help='The name of the effect/test allele',required='True')
    parser.add_argument('-a2', '--a2', help='The name of the REF or A2 (non-test allele)',required='True')
    parser.add_argument('-snp', '--snp', help='The name of the RSID field - this is required',required='True')
    parser.add_argument('-b', '--beta', help='The name of the beta column - if you have ORs, set this to -9',required='True')
    parser.add_argument('-OR', '--odds', help='If you are using ORs put it name of the column here',default=None)
    parser.add_argument('-p', '--pval', help='The name of the P-value column',required='True')

    results = parser.parse_args(args)
    return (results.input , results.output , results.a1 , results.a2 , results.snp , results.beta , results.odds, results.pval)


def main(input_file , output_file , a1 , a2 , snp , beta , OR, P): 
    raw=pd.read_csv(input_file, sep="\t")
    if beta == -9 : 
        temp = raw[[snp, a1, a2 , OR, P ]].copy()
        temp["BETA"] = np.log(temp[OR])
        temp = temp.drop(OR)
        temp.columns = ["SNP", "A1", "A2", "P", "BETA"]
        temp=temp.dropna()
        temp.to_csv(output_file, sep="\t", index=None)
    else: 
        temp = raw[[snp, a1, a2, beta, P]].copy()
        temp.columns=["SNP", "A1", "A2", "BETA", "P"]
        temp=temp.dropna()
        temp.to_csv(output_file,sep="\t", index=None)


if __name__ == '__main__':
  input_file , output_file , a1 , a2 , snp , beta , OR, P = check_arg(sys.argv[1:])
  main(input_file , output_file , a1 , a2 , snp , beta , OR, P )