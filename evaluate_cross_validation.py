# @ Project: TIMokokok
# @ Author: Mang Liang
# @ Date: 2024/11/14
# @ File content: describe the logic of the code

import os
import pandas as pd
from sklearn.model_selection import RepeatedKFold
import time
import warnings
import argparse

warnings.filterwarnings('ignore')

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



tpm = pd.read_csv('{}/{}'.format(work_path, ge), sep=' ')
rfk = RepeatedKFold(n_splits=5, n_repeats=1, random_state=87)
k = 0
for train_index, test_index in rfk.split(tpm):
    tpm_train = tpm.iloc[train_index, :]
    tpm_train.to_csv('{}/tpm_train.txt'.format(work_path), sep=' ', index=None)
    ge_cv = 'tpm_train.txt'
    t = time.localtime()
    print('locate time: {}-{}-{}, {}:{}:{}'.format(t[0], t[1], t[2], t[3], t[4], t[5]), )
    os.system('python genoSplit.py --work_path {} --geno {} --gtf {} --ge {} --window {} --core {}'.
              format(work_path, geno, gtf, ge_cv, window, 38))
    chr_files = os.listdir(work_path+'/TMP/genoChr')
    for file in chr_files:
        os.system('python dataProcess.py --work_path {} --geno {} --chr {} --gtf {} --ge {} --window {} --core {}'.
                      format(work_path, geno, file, gtf, ge_cv, window, core))
    t = time.localtime()
    print('locate time: {}-{}-{}, {}:{}:{}'.format(t[0], t[1], t[2], t[3], t[4], t[5]), )
    os.system('python single.py --work_path {} --core {} --predict {}'.format(work_path, core, 'False'))
    os.system('echo 123 | sudo -S sysctl -w vm.drop_caches=3 ')
    t = time.localtime()
    print('locate time: {}-{}-{}, {}:{}:{}'.format(t[0], t[1], t[2], t[3], t[4], t[5]), )
    os.system('python corMatrix.py --work_path {} --core {}'.format(work_path, core, ))
    t = time.localtime()
    print('locate time: {}-{}-{}, {}:{}:{}'.format(t[0], t[1], t[2], t[3], t[4], t[5]), )
    os.system('python multi.py --work_path {} --core {} --predict {}'.format(work_path, core, 'False'))
    # os.system('echo 123 | sudo -S sysctl -w vm.drop_caches=3 ')










