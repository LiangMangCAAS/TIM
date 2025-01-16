# @ Project: TIM
# @ Author: Mang Liang
# @ Date: 2024/11/14
# @ File content: describe the logic of the code

import pandas as pd
import numpy as np
from sklearn.model_selection import RepeatedKFold
from sklearn.linear_model import Ridge
# from cuml.linear_model import Ridge as cuRidge
import statsmodels.api as sm
# import cupy as cp
import warnings
warnings.filterwarnings('ignore')

def cal_rsquare(y, y_pre):
    X = sm.add_constant(y_pre)
    model = sm.OLS(y, X).fit()
    p = model.f_pvalue
    if p<0.05:
        return model.rsquared_adj
    else:
        return 0

# 根据基因组注释信息提取每个基因的位点信息，读取文件类型.gtf
def read_gtf(gtf):
    with open(gtf) as gtf:
        gtf_content = gtf.readlines()
    gtf.close()
    gene_location = []
    for i in range(5, len(gtf_content)):
        line_i = gtf_content[i].split()
        gene_location.append(line_i)
    gene_location = pd.DataFrame(gene_location)
    gene_location = gene_location.iloc[:, [0, 2, 3, 4, 8, 9]]
    gene_location = gene_location.loc[gene_location.iloc[:, 1] == 'gene']
    gene_location = gene_location.iloc[:, [0, 2, 3, 5]]
    gene_location.columns = ['chr', 'start', 'end', 'gene_id']
    gene_location['gene_id'] = gene_location['gene_id'].apply(lambda x: x[1:-2])
    return gene_location

def read_gff(gtf):
    with open(gtf) as gff:
        gff_content = gff.readlines()
    gff.close()
    gene_location = []
    for i in range( len(gff_content)):
        line_i = gff_content[i].split()
        if '#' not in line_i:
            gene_location.append(line_i)
    gene_location = pd.DataFrame(gene_location)
    gene_location = gene_location.iloc[:, [0, 2, 3, 4, 8, 9]]
    gene_location = gene_location.loc[gene_location.iloc[:, 1] == 'gene']
    gene_location_name = gene_location.iloc[:, 4].str.split(';',expand=True).iloc[:, 0].apply(lambda x:x[8:])
    gene_location = gene_location.iloc[:, [0, 2, 3, 4]]
    gene_location.columns = ['chr', 'start', 'end', 'gene_id']
    gene_location.iloc[:, -1] = gene_location_name
    return gene_location

def change_raw(data_name):
    with open(data_name) as f:
        file = f.readlines()
    res = []
    k= 0
    for i in file:
        res.append(i.split(' '))
        k+=1
    snp = pd.DataFrame(res)
    snp_id = snp.iloc[0, 6:]
    ind = snp.iloc[:, 1]
    snp = snp.iloc[1:, :]
    snp = snp.iloc[:,6:]
    snp = snp.replace('NA\n', 0)
    snp = snp.replace('NA', 0)
    snp = snp.astype(int)
    snp.index = ind[1:]
    bim = pd.read_csv('{}.bim'.format(data_name[:-4]), sep='\t', header=None)
    snp.columns = bim.iloc[:, 1]
    return snp


def tune_hp(x, y, cross_folds=5, repeat_num=1, random_seed=0,factor=1000, ):
    rfk = RepeatedKFold(n_splits=cross_folds, n_repeats=repeat_num, random_state=random_seed)
    k = 0
    matrix_pre = pd.DataFrame(np.zeros([x.shape[0],20]))
    alphas = np.array([((np.e / 2) ** i) / factor for i in np.arange(1, 41, 4)])
    accs = np.zeros([len(alphas), cross_folds * repeat_num])
    for alpha in alphas:
        acc = []
        pre_y = np.zeros(len(y))
        for train_index, test_index in rfk.split(x, y):
            estimator = Ridge(alpha=alpha)
            clf = estimator
            clf.fit(x[train_index, :], y[train_index])
            pre = clf.predict(x[test_index, :])
            acc_lm = cal_rsquare( y[test_index], pre)
            acc.append(acc_lm)
            pre_y[test_index] = pre
        accs[k, :] = acc
        matrix_pre.iloc[:, k] = pre_y
        k += 1
    accs_mean = accs.mean(axis=1)
    best_alpha_index = np.argsort(accs_mean)[-1]
    best_alpha = alphas[best_alpha_index]
    best_acc = accs_mean.max()
    best_pre = matrix_pre.iloc[:, best_alpha_index]
    return best_acc, best_alpha, best_pre.values



def interact_gene_pre(assist_gene, target_gene, ger, predict='False', factor=1000):
    matrix_pre = pd.DataFrame(np.zeros([ger.shape[0], len(assist_gene)]))
    matrix_pre.columns = assist_gene
    matrix_pre.index = ger.index
    x_target = pd.read_pickle('TMP/kernel_matrix/{}'.format(target_gene))
    vaild_matrix_pre = pd.DataFrame(np.zeros([x_target.shape[0] - ger.shape[0], len(assist_gene)]))
    # 利用五折交叉验证预测出所有个体的基因表达量预测值
    if predict=='True':
        vaild_matrix_pre = pd.DataFrame(np.zeros([x_target.shape[0]-ger.shape[0], len(assist_gene)]))
        x_target = x_target.drop(ger.index, axis=0)
        vaild_matrix_pre.index = x_target.index
        vaild_matrix_pre.columns = assist_gene
        for g in assist_gene:
            x_all =pd.read_pickle('TMP/kernel_matrix/{}'.format(g))
            x = x_all.loc[ger.index].values
            y = ger[target_gene].values
            acc, best_alpha, best_pre = tune_hp(x, y)
            matrix_pre[g] = best_pre
            rr = Ridge(alpha=best_alpha)
            rr.fit(x, y)
            vaild_pre = rr.predict(x_all.drop(ger.index, axis=0))
            vaild_matrix_pre[g] = vaild_pre
        return matrix_pre, vaild_matrix_pre
    else:
        for g in assist_gene:
            x_all =pd.read_pickle('TMP/kernel_matrix/{}'.format(g))
            x = x_all.loc[ger.index].values
            y = ger[target_gene].values
            acc, best_alpha, best_pre = tune_hp(x, y)
            matrix_pre[g] = best_pre
        return matrix_pre, vaild_matrix_pre








