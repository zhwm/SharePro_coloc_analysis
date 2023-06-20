import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='Prepare for SharePro simulation')
parser.add_argument('--dir', type=str, default=None, help='path to fastGWA result')
parser.add_argument('--loci', type=str, default=None, help='name of loci')
parser.add_argument('--ite', type=int, default=50, help='number of replications')
args = parser.parse_args()

for i in range(1, args.ite + 1):
    Qss = pd.read_csv(os.path.join(args.dir, 'Q{}.fastGWA'.format(i)), sep='\s+')
    Gss = pd.read_csv(os.path.join(args.dir, 'G{}.fastGWA'.format(i)), sep='\s+')
    Qss['Z'] = Qss['BETA'] / Qss['SE']
    Gss['Z'] = Gss['BETA'] / Gss['SE']
    Qss[['SNP', 'Z']].to_csv(os.path.join(args.dir, 'Q{}.z'.format(i)), sep='\t', index=False, header=False)
    Gss[['SNP', 'Z']].to_csv(os.path.join(args.dir, 'G{}.z'.format(i)), sep='\t', index=False, header=False)
    Qss[['SNP', 'A1', 'A2', 'AF1', 'BETA', 'SE', 'P', 'N']].to_csv(os.path.join(args.dir, 'Q{}.bse'.format(i)),
                                                                   sep='\t', index=False, header=True)
    Gss[['SNP', 'A1', 'A2', 'AF1', 'BETA', 'SE', 'P', 'N']].to_csv(os.path.join(args.dir, 'G{}.bse'.format(i)),
                                                                   sep='\t', index=False, header=True)
SHzld = pd.DataFrame({'z': ['Q{}.z,G{}.z'.format(i, i) for i in range(1, args.ite + 1)]})  # zld file for sharepro
SHzld['ld'] = '../{}.ld,../{}.ld'.format(args.loci, args.loci)
SHzld.to_csv(os.path.join(args.dir, 'SH.zld'), sep='\t', index=False, header=True)
