# coding = utf-8
# author: QiChen
# version: v1.5.0
# modification date: 2020/7/2

import sys, os, shutil, argparse, csv
if hasattr(sys, 'frozen'):
    os.environ['PATH'] = sys._MEIPASS + ";" + os.environ['PATH']
from collections import Counter
import platform
import multiprocessing
from lib import kmc_read, kmer_matrix, kmer_features, sequence_assembly
from lib import projectlist_file as plf

system_platform = platform.system()

def get_parameters():
    usage_info = "KmerGO [optional options] -i <input_files_folder> -t <input_trait_information>"
    parser = argparse.ArgumentParser(prog='KmerGO', usage=usage_info)
    parser.add_argument('-i', dest='sample_path', required=True, type=str, help='sample files path')
    parser.add_argument('-t', dest='trait_information_path', required=True, type=str, help='a csv file path of trait information')
    parser.add_argument('-m', '--mode', dest='mode', type=int, default=0, help='mode: 0-catagorical, 1-continuous (default: 0)')
    parser.add_argument('-k', '--kmerlength', dest='k_value', type=int, default=40, help='k-mer length (k from 14 to 256; default: 40)')
    parser.add_argument('-ci', dest='ci_value', type=int, default=2, help='minimal K-mer occurring times (default: 2)')
    parser.add_argument('-cs', dest='cs_value', type=int, default=65535, help='maximal K-mer occurring times (default: 65535)')
    parser.add_argument('-n', dest='process_number', type=int, default=2, help='number of processes (default: 24)')
    parser.add_argument('-assl', dest='ass_l', type=float, default=0.8, help='when mode = 0, logical features ASS value (default: 0.8)')
    parser.add_argument('-p', dest='p_value', type=float, default=0.01, help='numeric(mode=0) or logical(mode=1) features rank sum test p threshold value (default: 0.01)')
    parser.add_argument('-assn', dest='ass_n', type=float, default=0.8, help='when mode = 0, numeric features logistic regression ASS value (default: 0.8)')
    parser.add_argument('-corr', dest='corr_value', type=float, default=0.8, help='when mode = 1, numeric features coefficient of association ρ threshold value (default: 0.8)')
    return parser.parse_args()

def Check_csv_validity(param):
    global TI_dic, GroupA_Name, GroupA_Number, GroupB_Name, GroupB_Number, Group_Number
    TI_dic = {}
    GroupA_Name = None
    GroupB_Name = None
    GroupA_Number = 0
    GroupB_Number = 0
    Group_Number = 0
    try:
        csv_f = open(param.trait_information_path, 'r')
        reader = csv.reader(csv_f)
        for row in reader:
            if row[1] == 'trait':
                continue
            else:
                TI_dic[row[0]] = row[1]
    except:
        print('Error! Can not read the CSV file!', flush=True)
        return -1
    if len(TI_dic) == 0:
        print('Error! No samples!', flush=True)
        return -1
    if param.mode == 0:
        counter_dic = Counter(TI_dic.values())
        if len(counter_dic) != 2:
            print('Error! Group number must be 2!', flush=True)
            return -1
        else:
            GroupA_Name = list(counter_dic.keys())[0]
            GroupA_Number = list(counter_dic.values())[0]
            GroupB_Name = list(counter_dic.keys())[1]
            GroupB_Number = list(counter_dic.values())[1]
            Group_Number = GroupA_Number + GroupB_Number
            return 0
    else:
        for v in TI_dic.values():
            try:
                float(v)
            except:
                print('Error! Trait must be a float type!', flush=True)
                return -1
        return 0

def KMC_GO(param):
    kmc_thread = kmc_read.KMC_Thread((param.k_value, param.ci_value,
                                      param.cs_value, param.sample_path,
                                      'kmer_countings'))
    kmc_thread.start()
    last_info = ''
    print('Step1: k-mer counting', flush=True)
    while True:
        if kmc_thread.status <= 0:
            if kmc_thread.status == 0:
                return 0
            else:
                if kmc_thread.loginfo != last_info:
                    print(kmc_thread.loginfo, flush=True)
                return -1
        else:
            if kmc_thread.loginfo != last_info:
                last_info = kmc_thread.loginfo
                print(last_info, flush=True)

def GM_GO(param):
    tiplen = 0
    if param.mode == 0:
        gm_thread = kmer_matrix.GM_Thread(('kmer_countings', 'kmer_matrix',
                                           param.process_number, GroupA_Name,
                                           GroupA_Number, TI_dic))
    else:
        gm_thread = kmer_matrix.GM_Thread(('kmer_countings', 'kmer_matrix',
                                           param.process_number, '', len(TI_dic), TI_dic))
    gm_thread.start()
    last_info = ''
    print('Step2: k-mer union', flush=True)
    while True:
        if gm_thread.status <= 0:
            try:
                shutil.rmtree('temp')
                for job in gm_thread.jobs:
                    job.terminate()
            except:
                pass
            if gm_thread.status == 0:
                print('\r', end='', flush=True)
                for i in range(tiplen):
                    print(' ', end='', flush=True)
                print('\rTotal:100%%', flush=True)
                return 0
            else:
                print(gm_thread.loginfo, flush=True)
                return -1
        else:
            if gm_thread.status == 2:
                error_status = gm_thread.detective_error()
                if error_status == 0:
                    temp_progress = gm_thread.detective_progress()
                    tipstr = 'Total: %.2f%%|' % (100 * sum(temp_progress) /
                                                  sum(gm_thread.filesize))
                    for k in range(gm_thread.process_number):
                        if gm_thread.block_size[k] == 0:
                            tipstr += 'P' + str(k + 1) + ': 100%|'
                        else:
                            tipstr += 'P' + str(k + 1) + ': %.2f%%|' % \
                                      (100 * temp_progress[k] / gm_thread.block_size[k])
                    if last_info != tipstr:
                        last_info = tipstr
                        print('\r', end='', flush=True)
                        for i in range(tiplen):
                            print(' ', end='', flush=True)
                        print('\r' + tipstr[:-1], end='', flush=True)
                        tiplen = len(tipstr[:-1])
                else:
                    if error_status == -1:
                        gm_thread.loginfo = 'Error! Missing some KMC result files.'
                    elif error_status == -2:
                        gm_thread.loginfo = 'Error! Missing the result folder.'
                    elif error_status == -3:
                        gm_thread.loginfo = 'Error! Something is wrong when files\' pointer move.'
                    elif error_status == -10:
                        gm_thread.loginfo = 'Parameters error? Processes break down.'
                    gm_thread.status = -9

def GF_GO(param):
    tiplen = 0
    if param.mode == 0:
        catagorical_mode = True
    else:
        catagorical_mode = False
    gf_thread = kmer_features.GF_Thread(('kmer_matrix', 'kmer_features',
                                         param.ass_l, param.p_value,
                                         param.ass_n, GroupA_Number,
                                         GroupB_Number, GroupA_Name,
                                         GroupB_Name, TI_dic,
                                         param.corr_value, catagorical_mode))
    gf_thread.start()
    last_info = ''
    print('Step3: k-mer filtering', flush=True)
    while True:
        if gf_thread.status <= 0:
            try:
                shutil.rmtree('temp')
                for job in gf_thread.jobs:
                    job.terminate()
            except:
                pass
            if gf_thread.status == 0:
                print('\r', end='', flush=True)
                for i in range(tiplen):
                    print(' ', end='', flush=True)
                print('\rTotal:100%%')
                return 0
            else:
                print(gf_thread.loginfo, flush=True)
                return -1
        else:
            if gf_thread.status == 2:
                error_status = gf_thread.detective_error()
                if error_status == 0:
                    temp_progress = gf_thread.detective_progress()
                    tipstr = 'Total: %.2f%%|' % (100 * sum(temp_progress) /
                                                  sum(gf_thread.filesize))
                    for k in range(gf_thread.files_number):
                        tipstr += 'P' + str(k + 1) + ': %.2f%%|' % \
                                  (100 * temp_progress[k] / gf_thread.filesize[k])
                    if last_info != tipstr:
                        last_info = tipstr
                        print('\r', end='', flush=True)
                        for i in range(tiplen):
                            print(' ', end='', flush=True)
                        print('\r' + tipstr[:-1], end='', flush=True)
                        tiplen = len(tipstr[:-1])
                else:
                    if error_status == -1:
                        gf_thread.loginfo = 'Missing matrix files.'
                    elif error_status == -2:
                        gf_thread.loginfo = 'Missing the result folder.'
                    gf_thread.status = -9

def KA_GO(param):
    ka_thread = sequence_assembly.KA_Thread(('kmer_features', 'contig_result'))
    ka_thread.start()
    print('Step4: k-mer assembly', flush=True)
    while True:
        if ka_thread.status <= 0:
            if ka_thread.status == 0:
                print('Done! Result files are stored in \"contig_result\".', flush=True)
                return 0
            else:
                print(ka_thread.loginfo, flush=True)
                return -1

if __name__ == "__main__":
    multiprocessing.freeze_support()
    param = get_parameters()
    try:
        os.mkdir('kmer_countings')
    except:
        pass
    try:
        os.mkdir('kmer_matrix')
    except:
        pass
    try:
        os.mkdir('kmer_features')
    except:
        pass
    try:
        os.mkdir('contig_result')
    except:
        pass
    r = Check_csv_validity(param)
    if r == 0:
        r = KMC_GO(param)
    if r == 0:
        r = GM_GO(param)
    if r == 0:
        r = GF_GO(param)
    if r == 0:
        KA_GO(param)
    