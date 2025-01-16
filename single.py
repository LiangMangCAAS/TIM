# @ Project: TIM
# @ Author: Mang Liang
# @ Date: 2024/11/14
# @ File content: describe the logic of the code
from functions import tune_hp
import time
import pandas as pd
import numpy as np
import os
from multiprocessing import Pool  # 导入多进程中的进程池
import argparse
import warnings
import time
from sklearn.linear_model import Ridge
from joblib import Parallel, delayed
from functions import *
warnings.filterwarnings('ignore')

# set the work path
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument("--work_path", type=str, default="0", )
parser.add_argument("--core", type=int, default=5)
parser.add_argument("--predict", type=str, default="False")



args = parser.parse_args()
work_path = args.work_path
core = args.core
predict = args.predict

os.chdir(work_path)
# load the gene expression of reference population
ger = pd.read_pickle('ger')
# load the gene name
gene_names = os.listdir('TMP/kernel_matrix')


def single_fun(i, predict):
    gene = gene_names[i]
    if i%1000==0:
        print(i)
    gene_kernel = pd.read_pickle('TMP/kernel_matrix/{}'.format(gene))
    x = np.array(gene_kernel.loc[ger.index])
    y = np.array(ger[gene])
    nona_index = np.where(y > -99999)[0]
    x = x[nona_index, :]
    y = y[nona_index]
    acc, best_alpha,pre = tune_hp(x, y)
    if predict=='True':
        rr = Ridge(alpha=best_alpha)
        rr.fit(x, y)
        x_pre = gene_kernel.drop(ger.index, axis=0)
        ge_pre = rr.predict(x_pre)
        ge_pre = pd.DataFrame(ge_pre)
        ge_pre.index = x_pre.index
        t = time.localtime()
        t = time.asctime(t)
        ge_pre.to_pickle('TMP/preMatrix_single/{}_{}'.format(gene, t))
    return acc, best_alpha, gene


if __name__ == '__main__':
    time1 = time.time()
    p = Pool(core)
    results = []
    for i in range(len(gene_names)):
        r = p.apply_async(single_fun, args=(i,predict))
        results.append(r)
    p.close()
    p.join()
    res = []
    for i in results:
        res.append(i.get())
    pd.DataFrame(res).to_pickle('TMP/gene_params')
    time2 = time.time()
    print('single: time consuming ',int(time2 - time1), 's.')








