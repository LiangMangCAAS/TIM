# @ Project: TIM
# @ Author: Mang Liang
# @ Date: 2024/11/14
# @ File content: describe the logic of the code

import pandas as pd
import os
import numpy as np
import argparse
from scipy.stats import pearsonr
from multiprocessing import Pool
import time
import warnings
# import cupy as cp
warnings.filterwarnings('ignore')

############################################
# 直接调用最佳的参数，进行每个基因表达量的预测
############################################

# load the gene_names

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument("--work_path", type=str, default="0", )
parser.add_argument("--core", type=int, default=5)
args = parser.parse_args()

work_path = args.work_path
core = args.core


os.chdir(work_path)
gene_names = os.listdir('TMP/kernel_matrix')
ger = pd.read_pickle('ger')
ger_matrix = ger[gene_names].values
assist_genes_all = gene_names

mean = np.mean(ger_matrix, axis=0)
std_dev = np.std(ger_matrix, axis=0)
# 中心化数据（减去均值）
matrix_centered = ger_matrix - mean
correlation_matrix = np.dot(matrix_centered.T, matrix_centered) / (ger_matrix.shape[0] - 1)
correlation_matrix = correlation_matrix / np.outer(std_dev, std_dev)
correlation_matrix = pd.DataFrame(correlation_matrix)
correlation_matrix.columns = assist_genes_all
correlation_matrix.index = assist_genes_all
t = time.localtime()
t = time.asctime(t)
correlation_matrix.to_pickle('TMP/correlation_matrix_{}'.format(t))


# 直接统一计算
def calcaulate_cor_matrix(i):
    gene = gene_names[i]
    correlation_matrix.loc[gene].to_pickle('TMP/corMatrix/{}'.format(gene))
    return gene


if __name__ == '__main__':
    time1 = time.time()
    p = Pool(core)
    results = []
    for i in range(len(gene_names)):
    # for i in range(100):
        r = p.apply_async(calcaulate_cor_matrix, args=(i,))
        results.append(r)
    p.close()
    p.join()
    time2 = time.time()
    print('corMatrix: time consuming ',int(time2 - time1), 's.')

