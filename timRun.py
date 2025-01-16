# @ Project: TIM
# @ Author: Mang Liang
# @ Date: 2024/11/26
# @ File content: describe the logic of the code

import os
import pandas as pd
import time
import warnings
warnings.filterwarnings('ignore')

work_path = 'example'
geno = 'genotype'
gtf = 'gene_annotation.gtf'
window = 1000*1000
n_core = 20
ge = 'adjusted_tpm.txt'
phenotype = 'phenotype.txt'
trait_name = 'trait_name'


t = time.localtime()
print('locate time: {}-{}-{}, {}:{}:{}'.format(t[0], t[1], t[2], t[3], t[4], t[5]), )
os.system('python evaluate_cross_validation.py --work_path {} --geno {} --gtf {} --ge {} --window {} --core {}'.format(work_path, geno, gtf, ge, window, n_core))

t = time.localtime()
print('locate time: {}-{}-{}, {}:{}:{}'.format(t[0], t[1], t[2], t[3], t[4], t[5]), )
os.system('python verify.py --work_path {} --ge {}'.format(work_path, ge))


t = time.localtime()
print('locate time: {}-{}-{}, {}:{}:{}'.format(t[0], t[1], t[2], t[3], t[4], t[5]), )
os.system('python association.py --work_path {} --trait_name {} --trait_name {}'.format(work_path, trait_name, phenotype ))




