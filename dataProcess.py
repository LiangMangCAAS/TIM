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
parser.add_argument("--chr", type=str, default='0')
parser.add_argument("--ge", type=str, default='0')
parser.add_argument("--window", type=int, default=500000)
parser.add_argument("--core", type=int, default=5)

# 提取参数
args = parser.parse_args()
work_path = args.work_path
gtf = args.gtf
geno = args.geno
chr_file = args.chr
ge = args.ge
window = args.window
core = args.core


# 检查设置的工作路径是否存在
if not os.path.exists(work_path):
    print('Error: Please check the work path !')
    os._exit(0)
# 切换工作路径

os.chdir(work_path)

if not os.path.exists(gtf):
    print('Error: Please check the file path of the gene annotation file !')
    os._exit(0)
if not os.path.exists(ge):
    print('Error: Please check the file path of the gene expression file !')
    os._exit(0)
if not os.path.exists('{}.fam'.format(geno)):
    print('Error: Please check the file path of the genome file !')
    os._exit(0)

# 读取注释文件
if 'gtf' in gtf:
    if 'gene_location' in os.listdir():
        gene_location = pd.read_pickle('gene_location')
    else:
        gene_location = read_gtf(gtf)
        gene_location.to_pickle('gene_location')
else:
    if 'gene_location' in os.listdir():
        gene_location = pd.read_pickle('gene_location')
    else:
        gene_location = read_gff(gtf)
        gene_location.to_pickle('gene_location')


# 读取转录组数据
ger = pd.read_csv(ge, sep=' ', index_col=0)


# 创建SNP位点信息矩阵
bim = pd.read_csv('{}.bim'.format(geno), sep='\t', header=None)
genotype_bim_df = pd.DataFrame(np.zeros([bim.shape[0], 3]))
genotype_bim_df.columns = ['snp_name', 'chr', 'pos']
genotype_bim_df['snp_name'] =bim.iloc[:, 1]
genotype_bim_df['chr'] = bim.iloc[:, 0].astype(str)
genotype_bim_df['pos'] =bim.iloc[:, 3]
bim.index = bim.iloc[:, 1]
chrs = np.unique(genotype_bim_df['chr'])
# genotype.columns = genotype_bim_df['snp_name'].values


num = []
for i in chrs:
    if len(i)>2:
        num.append(i)

chrs = np.delete(chrs, num)
chrs = chrs.astype(str)



genotype = pd.read_pickle('TMP/genoChr/{}'.format(chr_file ))


samples = genotype.index
sample_reference = np.intersect1d(samples, ger.index)
sample_predict = np.setdiff1d(samples, sample_reference)

id_resort = pd.DataFrame(np.zeros([len(sample_reference)+len(sample_predict), 1]))
id_resort.iloc[:len(sample_reference),0] = sample_reference
id_resort.iloc[len(sample_reference):,0] = sample_predict
id_resort.index = id_resort.iloc[:,0]

ger = ger.loc[sample_reference]
ger.to_pickle('ger')


# 创建函数获取每个基因的cissnp集合
def generate_gene_cissnp_kernel_matrix(i):
    gene_name = ger.columns[i]
    gene_information =  gene_location.loc[gene_location['gene_id']==gene_name]
    chr = gene_information['chr'].values[0]
    if ('chr'+str(chr)) != chr_file:
        return False
    else:
        start = int(gene_information['start'])-window
        end = int(gene_information['end']) + window
        if str(chr) not in chrs:
            return False
        else:
            bim_chr = genotype_bim_df.loc[genotype_bim_df['chr']==str(chr)]
            snp_select = bim_chr['pos'].apply(lambda x: start<int(x) & int(x)<end)
            if snp_select.sum() == 0:
                return False
            else:
                gene_cissnp_group = bim_chr.loc[snp_select]['snp_name'].values
                cissnp = genotype[gene_cissnp_group].values
                cissnp = pd.DataFrame(cissnp)
                cissnp.columns = gene_cissnp_group
                cissnp.index = samples
                cissnp = cissnp.fillna(0)
                cissnp = cissnp.loc[id_resort.index]
                if cissnp.shape[1] > 100000:
                    # print(gene_name)
                    return False
                kernel_cissnp = rbf_kernel(cissnp)[:, :len(ger.index)]
                kernel_cissnp = pd.DataFrame(kernel_cissnp)
                kernel_cissnp.index = id_resort.index
                kernel_cissnp.columns = id_resort.index[:len(ger.index)]
                kernel_cissnp.to_pickle('{}/TMP/kernel_matrix/{}'.format(work_path, gene_name))
                return gene_name


if __name__ == '__main__':
    time1 = time.time()
    p=Pool(core) #创建含有十个10个进程的进程池
    results=[] #存放每一个进程返回的结果
    for i in range(ger.shape[1]): # 启动10个进程
        r=p.apply_async(generate_gene_cissnp_kernel_matrix,args=(i,)) # 产生一个非同步进程，函数newsin的参数用args传递
        results.append(r) # 将返回结果放入results
    p.close() #关闭进程池
    p.join()  #结束
    time2 = time.time()


