import os
import statsmodels.stats.multitest as smm
from scipy.stats import chi2
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np
from statsmodels.api import OLS, add_constant
from sklearn.decomposition import PCA
import argparse
import warnings
warnings.filterwarnings('ignore')


parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument("--work_path", type=str, default="0", )
parser.add_argument("--trait_name", type=str, default="0", )
parser.add_argument("--phenotype", type=str, default="0", )

# 提取参数
args = parser.parse_args()
work_path = args.work_path
trait_name = args.trait_name
phenotype = args.phenotype

os.chdir(work_path)
phes = pd.read_excel(phenotype, sep=' ', index_col=0)

tran_files = os.listdir('results/final_imputation')
data_0 = pd.read_pickle('results/final_imputation/{}'.format(tran_files[0]))
trans = pd.DataFrame(np.zeros([data_0.shape[0], len(tran_files)])).astype(float)
trans.columns = tran_files
trans.index = data_0.index

for f in range(len(tran_files)):
    if f %1000==0:
        print(f)
    file = tran_files[f]
    tran_file = pd.read_pickle('results/final_imputation/{}'.format(file))
    trans.iloc[:, f] = tran_file.values.flatten()

trans.columns = tran_files
# trans.to_pickle('trans')
# trans = pd.read_pickle('trans')

pre_acc = pd.read_excel('results/accuracy_results.xlsx', index_col=0)
trans = trans.iloc[:, (np.where(pre_acc['multi_acc']>0.1))[0]]

commen_id = np.intersect1d(trans.index, phes.index)
trans = trans.loc[commen_id]
phes = phes.loc[commen_id]


# 2. 定义 GWAS 分析函数
def gwas_analysis(data, phenotype_col, gene_cols, ):
    results = []
    y = data[phenotype_col]  # 提取表型数据
    k=0
    for gene in gene_cols:
        X = add_constant(data[gene])
        model = OLS(y, X).fit()  # 线性回归
        p_value = model.pvalues[gene]
        beta = model.params[gene]
        results.append({"GENE": gene, "Beta": beta, "P-value": p_value})
        k+=1
    results_df = pd.DataFrame(results)
    return results_df

# 3. 运行 GWAS 分析

data = pd.concat([trans, phes], axis=1)
results_df = gwas_analysis(data, trait_name, trans.columns,)

results_df['P_value_FDR'] = smm.multipletests(results_df['P-value'], method='fdr_bh')[1]
# 4. 标记显著性基因
results_df['Significant'] = results_df['P_value_FDR'] < 0.05
results_df.to_excel('results/TIM_results.xlsx')

