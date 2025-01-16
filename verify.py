import pandas as pd
import numpy as np
import os
import statsmodels.api as sm
import argparse
import warnings
from functions import cal_rsquare

warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument("--work_path", type=str, default="0", )
parser.add_argument("--ge", type=str, default='0')

args = parser.parse_args()
work_path = args.work_path
ge = args.ge

os.chdir(work_path)

def vertify_acc(workPath):
    files = os.listdir(workPath)
    files_index = np.array([i[:18] for i in files])
    genes = np.unique(files_index)
    files = pd.DataFrame(files)
    files.index = files_index
    acc_all = []
    tpm_true = pd.read_csv(ge, sep=' ', index_col=0)
    for gene in genes:
        gene_file = files.loc[gene].values.flatten()
        accs = []
        for i in gene_file:
            gene_pre = pd.read_pickle('{}/{}'.format(workPath, i))
            gene_true = tpm_true.loc[gene_pre.index][gene]
            acc = cal_rsquare(gene_true.values, gene_pre.values)
            accs.append(acc)
        gene_acc = np.array(accs).mean()
        acc_all.append(gene_acc)
    acc_all = pd.DataFrame(acc_all)
    acc_all.index = genes
    return acc_all


os.chdir(work_path)
tpm_true = pd.read_csv(ge, sep=' ', index_col=0)
single_acc = vertify_acc('TMP/preMatrix_single')
multi_acc = vertify_acc('TMP/preMatrix_multi')
print( '   single_acc: ', single_acc.mean().values)
print( '   multi_acc: ', multi_acc.mean().values)
files = os.listdir('TMP/preMatrix_multi')
files_index = np.array([i[:18] for i in files])
genes = np.unique(files_index)
cors = []
accs =[]
for i in genes:
    cor = pd.read_pickle('TMP/corMatrix/{}'.format(i)).sort_values(ascending=False)
    cor_mean = cor[:int(cor.shape[0]*0.01)].mean()
    acc_mean = single_acc.loc[cor[:int(cor.shape[0]*0.01)].index].mean().values[0]
    cors.append(cor_mean)
    accs.append(acc_mean)

cors = pd.DataFrame(cors)
cors.index = genes

accs = pd.DataFrame(accs)
accs.index = genes

commen_index = np.intersect1d(single_acc.index, multi_acc.index)

results = pd.DataFrame(np.zeros([len(commen_index), 4]))
results.index = commen_index
results.columns = ['cor_mean', 'single_acc', 'multi_acc', 'cor_acc_mean']
results.iloc[:, 0] = cors.loc[commen_index]
results.iloc[:, 1] = single_acc.loc[commen_index]
results.iloc[:, 2] = multi_acc.loc[commen_index]
results.iloc[:, 3] = accs.loc[commen_index]

results.to_excel('results/accuracy_results.xlsx')












