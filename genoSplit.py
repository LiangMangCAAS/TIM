# @ Project: TIM
# @ Author: Mang Liang
# @ Date: 2024/11/14
# @ File content: describe the logic of the code


import argparse
import os
from multiprocessing import Pool
from sklearn.metrics.pairwise import rbf_kernel
import time
from functions import *
import warnings
import pandas as pd
import numpy as np

warnings.filterwarnings('ignore')

#
#
# 设置代码接收的参数
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument("--work_path", type=str, default="0", )
parser.add_argument("--gtf", type=str, default="0", )
parser.add_argument("--geno", type=str, default='0')
parser.add_argument("--ge", type=str, default='0')
parser.add_argument("--window", type=int, default=500000)
parser.add_argument("--core", type=int, default=5)

# 提取参数
args = parser.parse_args()
work_path = args.work_path
gtf = args.gtf
geno = args.geno
ge = args.ge
window = args.window
core = args.core

# 检查设置的工作路径是否存在
if not os.path.exists(work_path):
    print('Error: Please check the work path !')
    os._exit(0)
# 切换工作路径

os.chdir(work_path)

# 创建临时文件夹和结果文件夹
if not os.path.exists('TMP'):
    os.makedirs('TMP/kernel_matrix')
    os.makedirs('TMP/corMatrix')
    os.makedirs('TMP/preMatrix_single')
    os.makedirs('TMP/preMatrix_multi')
    os.makedirs('TMP/genoChr')

if not os.path.exists('results'):
    os.makedirs('results/final_imputation')
    os.makedirs('results/GWAS_results')



# 读取转录组数据
ger = pd.read_csv(ge, sep=' ', index_col=0)

# 创建SNP位点信息矩阵
bim = pd.read_csv('{}.bim'.format(geno), sep='\t', header=None)
genotype_bim_df = pd.DataFrame(np.zeros([bim.shape[0], 3]))
genotype_bim_df.columns = ['snp_name', 'chr', 'pos']
genotype_bim_df['snp_name'] = bim.iloc[:, 1]
genotype_bim_df['chr'] = bim.iloc[:, 0].astype(str)
genotype_bim_df['pos'] = bim.iloc[:, 3]
bim.index = bim.iloc[:, 1]
chrs = np.unique(genotype_bim_df['chr'])
# genotype.columns = genotype_bim_df['snp_name'].values


num = []
for i in chrs:
    if len(i) > 2:
        num.append(i)

chrs = np.delete(chrs, num)
chrs = chrs.astype(str)

for i in chrs:
    print(i)
    os.system(
        '/media/mang/diskB/./plink --bfile {} --chr {} --make-bed --recode A --out TMP/chr{} --chr-set 50'.format(geno,
                                                                                                                  i, i))
    genotype_chr = change_raw('TMP/chr{}.raw'.format(i))
    genotype_chr.to_pickle('TMP/genoChr/chr{}'.format(i))



