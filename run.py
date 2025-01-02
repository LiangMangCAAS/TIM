# @ Project: TIM
# @ Author: Mang Liang
# @ Date: 2024/11/14
# @ File content: describe the logic of the code

import os
import pandas as pd
from sklearn.model_selection import RepeatedKFold
import time
import warnings

warnings.filterwarnings('ignore')
species = 'cattle335'
tissue = 'muscle'
print(species, tissue)
work_path = '/media/mang/diskA/geim/{}/{}'.format(species, tissue)
geno = 'genotype'
gtf = 'gene_annotation.gtf'
window = 500*1000 # 1MB
ge = 'adjusted_tpm.txt'


# '/media/mang/diskA/geim/{}/{}' 'C:/work/tim/{}'
for tissue in os.listdir('/media/mang/diskA/geim/{}/'.format(species, )):
    work_path = '/media/mang/diskA/geim/{}/{}'.format(species, tissue)
    # t = time.localtime()
    # print('locate time: {}-{}-{}, {}:{}:{}'.format(t[0], t[1], t[2], t[3], t[4], t[5]), )
    # os.system('python genoSplit.py --work_path {} --geno {} --gtf {} --ge {} --window {} --core {}'.
    #           format(work_path, geno, gtf, ge, window, 38))
    #
    # chr_files = os.listdir(work_path+'/TMP/genoChr')
    # for file in chr_files:
    #     os.system('python dataProcess.py --work_path {} --geno {} --chr {} --gtf {} --ge {} --window {} --core {}'.
    #                   format(work_path, geno, file, gtf, ge, window, 38))
    #
    # t = time.localtime()
    # print('locate time: {}-{}-{}, {}:{}:{}'.format(t[0], t[1], t[2], t[3], t[4], t[5]), )
    # os.system('python single.py --work_path {} --core {} --predict {}'.format(work_path, 38, True))
    #
    # t = time.localtime()
    # print('locate time: {}-{}-{}, {}:{}:{}'.format(t[0], t[1], t[2], t[3], t[4], t[5]), )
    # os.system('python corMatrix.py --work_path {} --core {}'.format(work_path, 78, ))

    t = time.localtime()
    print('locate time: {}-{}-{}, {}:{}:{}'.format(t[0], t[1], t[2], t[3], t[4], t[5]), )
    os.system('python multi.py --work_path {} --core {} --predict {}'.format(work_path, 20, True))












