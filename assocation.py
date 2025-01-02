import os
import statsmodels.stats.multitest as smm
from scipy.stats import chi2
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np
from statsmodels.api import OLS, add_constant
from sklearn.decomposition import PCA



work_path = 'C:\\work\\cattle'
os.chdir(work_path)
data_name = 'genotype/genotype.raw'

# load genotype
# snp = change_raw(data_name,)
# snp.to_pickle('genotype/genotype')
snp = pd.read_pickle('genotype/genotype')

# load transcriptome data
# tran_files = os.listdir('preMatrix')
#
# data_0 = pd.read_pickle('preMatrix/{}'.format(tran_files[0]))
# trans = pd.DataFrame(np.zeros([data_0.shape[0], len(tran_files)])).astype(float)
# trans.columns = tran_files
# trans.index = data_0.index
#
# for f in range(len(tran_files)):
#     if f %1000==0:
#         print(f)
#     file = tran_files[f]
#     tran_file = pd.read_pickle('preMatrix/{}'.format(file))
#     trans.iloc[:, f] = tran_file.values.flatten()
#
# trans.columns = tran_files
# trans.to_pickle('trans')
trans = pd.read_pickle('trans')

commen_individial = np.intersect1d(trans.index, snp.index)

trans = trans.loc[commen_individial]
snp = snp.loc[commen_individial]


pre_acc = pd.read_excel('accuracy_Sat Dec 14 02_55_38 2024.xlsx', index_col=0)
pre_acc['max_acc'] = pre_acc.max(axis=1)
trans = trans.iloc[:, (np.where(pre_acc['max_acc']>0.1))[0]]


# correct the phe
phe = pd.read_excel('屠宰数据-修正后（2024-1-23）.xlsx')
fix_effect = ['固定效应-屠宰日期', '固定效应-性别', '固定效应-育肥场', '育肥期天数（修正后）','固定效应-屠宰场', ]
target_trait = ['宰前活重', '育肥期日增重（修正后）', '胴体重（修正后）', '骨重（修正后）', '眼肌面积（12肋）']

phe['背最长肌'] = phe['外脊（后）']+phe[  '眼肌（后+取样）(修正后)']+phe[ '上脑（后）（修改后）']

phe.index = phe.iloc[:, 0]
phe = phe.loc[snp.index]

np.random.seed(87)
def correct_phe(phe, fix_effect, target_trait, tans):
    data_fix = phe[fix_effect]
    data_p = phe[target_trait].astype(float)
    data_p = data_p.dropna()
    data_fix = data_fix.dropna(how='any')
    commen_id = np.intersect1d(data_fix.index, data_p.index)
    data_p = data_p.loc[commen_id]
    data_fix = data_fix.loc[commen_id]
    snp_data = tans.loc[commen_id]
    pca = PCA(n_components=5)  # 提取前 10 个主成分
    principal_components = pca.fit_transform(snp_data)
    principal_components = pd.DataFrame(principal_components)
    principal_components.index = snp_data.index
    data_fix = pd.concat([data_fix, principal_components], axis=1)
    X = sm.add_constant(data_fix)  # 添加截距项
    y = data_p
    model = sm.OLS(y, X, ).fit()
    return y - model.predict(X)

fix_effect = ['固定效应-屠宰日期', '固定效应-性别',  '屠宰日龄(修正后)','固定效应-屠宰场','育肥期天数（修正后）', '胴体重（修正后）',]# '胴体重（修正后）'

target_trait =   '上脑（后）（修改后）' # 眼肌（后+取样）(修正后)     胴体重（修正后） 上脑（后）（修改后）   背最长肌  外脊（后）


lw_correct = correct_phe(phe, fix_effect, target_trait, trans)

commen_index = np.intersect1d(snp.index, lw_correct.index)
commen_index = np.intersect1d(commen_index, trans.index)


data_fix = phe[fix_effect].loc[lw_correct.index]
lw_correct = pd.DataFrame(lw_correct)
lw_correct.columns = [target_trait]
data = pd.concat([trans.loc[commen_index], lw_correct, data_fix], axis=1)



# 2. 定义 GWAS 分析函数
def gwas_analysis(data, phenotype_col, genotype_cols, ):
    results = []
    y = data[phenotype_col]  # 提取表型数据
    k=0
    # fix = principal_components
    # X_cov = add_constant(fix)  # 协变量设计矩阵
    for snp in genotype_cols:
        if k % 1000 == 0:
            print(k)
        # X = pd.concat([X_cov, data[snp]], axis=1)  # 添加 SNP 到设计矩阵
        X = add_constant(data[snp])
        model = OLS(y, X).fit()  # 线性回归
        p_value = model.pvalues[snp]
        beta = model.params[snp]
        results.append({"SNP": snp, "Beta": beta, "P-value": p_value})
        k+=1
    results_df = pd.DataFrame(results)
    return results_df


# 3. 运行 GWAS 分析
phenotype_col = target_trait
genotype_cols = trans.columns
results_df = gwas_analysis(data, phenotype_col, genotype_cols,)





results_df['P_value_FDR'] = smm.multipletests(results_df['P-value'], method='fdr_bh')[1]
# 4. 标记显著性基因
results_df['Significant'] = results_df['P_value_FDR'] < 0.05

# 输出结果
# print(results_df.loc[results_df['Significant']==True])
# print(target_trait)
# print(results_df['P_value_FDR'].describe())
# print(results_df['P-value'].describe())




gene_location = pd.read_pickle('gene_location')
gene_location.index = gene_location.iloc[:, -1]

results_df.index = results_df.iloc[:, 0]
results_df['chr'] = gene_location.loc[results_df.index]['chr']
results_df['pos'] = gene_location.loc[results_df.index]['end'].astype(int)-500000

import geneview as gv
import matplotlib.pyplot as plt
from geneview.gwas import manhattanplot,qqplot

ax = qqplot(results_df["P-value"], color="#00bb33", xlabel="Expected p-value(-log10)", ylabel="Observed p-value(-log10)")
plt.savefig('cattle_QQ_plot_{}'.format(target_trait),  dpi=300, bbox_inches='tight')


data = results_df[['chr', 'SNP', 'pos', 'P-value']]
data.columns = ['#CHROM','rsID','POS','P']
data['#CHROM'] = data['#CHROM'].astype(int)
data['POS'] = data['POS'].astype(int)
data = data.sort_values(by=['#CHROM', 'POS'])
ax = manhattanplot(data, xlabel="Chromosome", ylabel="-Log10(P-value)",
                   genomewideline=-np.log10(0.05/data.shape[0]),)  # 这就是Manhattan plot的函数
plt.savefig('manhattan_plot_{}'.format(target_trait),  dpi=300, bbox_inches='tight')

sig_snp = results_df.loc[results_df['P-value']<(0.00001)]
print(sig_snp)
print(results_df.loc[results_df['Significant']==True])
print(target_trait)
results_df["Chi2_Observed"] = chi2.ppf(1 - results_df["P-value"], df=1)
# Step 2: 计算膨胀系数 λGC
chi2_median_observed = np.median(results_df["Chi2_Observed"])
chi2_median_expected = chi2.ppf( 0.5, df=1)  # 自由度为 1，理论中位数 = 0.4549
lambda_gc = chi2_median_observed / chi2_median_expected
print('P_values膨胀系数： ', lambda_gc)


results_df.loc[results_df['Significant']==True].to_excel('significant_gene_{}.xlsx'.format(target_trait))

gene_location.loc['ENSBTAG00000002126']