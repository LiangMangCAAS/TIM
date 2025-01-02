# @ Project: TIM
# @ Author: Mang Liang
# @ Date: 2024/11/14
# @ File content:
from functions import *
import pandas as pd
import os
import numpy as np
import argparse
from multiprocessing import Pool
import time
import warnings
from sklearn.linear_model import Ridge
from functions import tune_hp, cal_rsquare
warnings.filterwarnings('ignore')
import time

############################################
# 直接调用最佳的参数，进行每个基因表达量的预测
############################################

# load the gene_names

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument("--work_path", type=str, default="0", )
parser.add_argument("--core", type=int, default=5)
parser.add_argument("--predict", type=bool, default=False)

args = parser.parse_args()
work_path = args.work_path
core = args.core
predict = args.predict

os.chdir(work_path)

gene_names = os.listdir('TMP/kernel_matrix')
ger = pd.read_pickle('ger')
ger = ger[gene_names]



# load the best parameters of the Ridge
gene_params = pd.read_pickle('TMP/gene_params')
gene_params.columns = ['acc', 'best_alpha', 'gene_id']
gene_params.index = gene_params.iloc[:,-1].values
# 直接统一计算


def multi_funB(i):
    gene = gene_names[i]
    y = ger[gene].values
    cors = pd.read_pickle('TMP/corMatrix/{}'.format(gene))
    assist_gene = cors.sort_values()[::-1][:int(len(cors)*0.01)].index
    x, x_pre = interact_gene_pre(assist_gene, gene, ger)
    acc, best_alpha, pre = tune_hp(x.values, y, factor=100)
    original_acc =  gene_params.iloc[:, 0].loc[gene]
    if (acc > original_acc) and predict:
        rr = Ridge(alpha=best_alpha)
        rr.fit(x.values, y)
        ge_pre = rr.predict(x_pre)
        ge_pre = pd.DataFrame(ge_pre)
        ge_pre.index = x_pre.index
        ge_pre.to_pickle('TMP/preMatrix/{}'.format(gene))
    else:
        acc = original_acc
    if i % 1000 == 0:
        t = time.localtime()
        print('####:   ', i, '    ####, ', 'locate time: {}-{}-{}, {}:{}:{}'.format(t[0], t[1], t[2], t[3], t[4], t[5]), )
        print('gene: ', gene, '   original_acc: ', gene_params.iloc[:, 0].loc[gene], '   planB_acc: ', acc)
    return gene, gene_params.iloc[:, 0].loc[gene], acc,




if __name__ == '__main__':
    time1 = time.time()
    p = Pool(core)
    results = []
    # for i in range(100):
    for i in range(len(gene_names)):
        r = p.apply_async(multi_funB, args=(i,))
        results.append(r)
    p.close()
    p.join()
    res1 = pd.DataFrame(np.zeros([len(gene_names), 2]))
    res1.columns = ['original_acc', "multi_acc"]
    res1.index = gene_names[:len(gene_names)]
    for i in results:
        res_i = i.get()
        res1.loc[res_i[0]] = [res_i[1], res_i[2]]
    time2 = time.time()
    t = time.localtime()
    t = time.asctime(t)
    res1.to_excel('results/accuracy_{}.xlsx'.format(t))
    print('multi: time consuming ',int(time2 - time1), 's.')




