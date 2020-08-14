# coding = utf-8
# author: QiChen
# version: v5.0
# modification date: 2020/7/23

import os
import time
import threading
import _thread
from scipy import stats
import numpy as np
import pandas as pd
from multiprocessing import Process
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix

def Categorical_feature_filtering(Nprocess, input_path, output_path1, output_path2, param):
    Number_of_A_Group = param[0]
    Number_of_B_Group = param[1]
    A_Name = param[2]
    B_Name = param[3]
    TI_dic = param[4]
    ass_l = param[5]
    wicxon_p = param[6]
    ass_n = param[7]
    bufsize = param[8]

    try:
        fi = open(input_path, 'rb')
    except:
        ferr = open(os.path.join('temp', 'GF_error_status=-1'), 'w')
        ferr.close()
        return
    try:
        fo1 = open(output_path1,'wb', buffering=bufsize)
        fo2 = open(output_path2, 'wb', buffering=bufsize)
    except:
        ferr = open(os.path.join('temp', 'GF_error_status=-2'), 'w')
        ferr.close()
        return

    # variate initialization
    last_progress = 0
    headtext = fi.readline().decode('utf-8')
    progress = len(headtext)
    headtext = headtext.strip()
    headlist = headtext.split('\t')
    headlist = headlist[1:]
    strline = fi.readline().decode('utf-8')

    # write head text
    fo1.write((headtext + '\tASS\tLabel\n').encode('utf-8'))
    fo2.write((headtext + '\tASS-l\tP\tASS-n\tLabel\n').encode('utf-8'))

    # start time
    last_time = time.time()

    # main loop
    while strline != '':
        st = strline.find('\t') + 1
        si = strline[st:-1]
        si = si.split('\t')
        tp = 0
        fn = 0
        group1 = []
        group2 = []
        label_Y = []
        # logical filtering
        for i in range(len(si)):
            if not(headlist[i] in TI_dic):
                continue
            nor_value = float(si[i])
            if TI_dic[headlist[i]] == A_Name:
                group1.append(nor_value)
                label_Y.append(0)
                if si[i] != '0':
                    tp += 1
            else:
                group2.append(nor_value)
                label_Y.append(1)
                if si[i] != '0':
                    fn += 1
        ASStrue = tp / Number_of_A_Group + (Number_of_B_Group - fn) / Number_of_B_Group
        ASStrue /= 2
        if ASStrue > 1 - ASStrue:
            label = A_Name
        else:
            ASStrue = 1 - ASStrue
            label = B_Name
        if ASStrue >= ass_l:
            fo1.write(strline[:-1].encode('utf-8'))
            fo1.write(('\t' + str(ASStrue) + '\t' + label + '\n').encode('utf-8'))
        else:
            # numerical filtering
            (Zvalue, Pvalue) = stats.ranksums(group1, group2)
            if Pvalue < wicxon_p:
                train_X = np.array(group1 + group2)
                train_X = train_X.reshape(len(train_X), 1)
                excepted = train_Y = np.array(label_Y)
                model = LogisticRegression(C=1e9)
                model.fit(train_X, train_Y)
                predicted = model.predict(train_X)
                ConfusionMtrix = confusion_matrix(excepted, predicted, labels=[0, 1])   # 0-healthy 1-patient
                cc = 0.5 * (float(ConfusionMtrix[0, 0]) / float(Number_of_A_Group) + float(ConfusionMtrix[1, 1]) /
                            float(Number_of_B_Group))
                if cc >= ass_n:
                    if np.mean(group1) > np.mean(group2):
                        label = A_Name
                    else:
                        label = B_Name
                    fo2.write(strline[:-1].encode('utf-8'))
                    fo2.write(('\t' + str(ASStrue) + '\t' + str(Pvalue) + '\t' + str(cc) + '\t' + label + '\n').encode('utf-8'))
        progress += len(strline)     # update the progress
        if time.time() - last_time >= 1:    # output the progress file per second
            os.rename(os.path.join('temp', 'GF_progress' + str(Nprocess) + ' ' + str(last_progress)),
                      os.path.join('temp', 'GF_progress' + str(Nprocess) + ' ' + str(progress)))
            last_progress = progress
            last_time = time.time()
        strline = fi.readline().decode('utf-8')

    # progress 100%
    os.rename(os.path.join('temp', 'GF_progress' + str(Nprocess) + ' ' + str(last_progress)),
              os.path.join('temp', 'GF_progress' + str(Nprocess) + ' ' + str(progress)))

    # close file
    fi.close()
    fo1.close()
    fo2.close()

    # create ok-flag file
    fok = open(os.path.join('temp', str(Nprocess) + '_ok'), 'w')
    fok.close()

def Continuous_feature_filtering(Nprocess, input_path, output_path1, output_path2, param):
    TI_dic = param[0]
    wicxon_p = param[1]
    corr_value = param[2]
    bufsize = param[3]

    try:
        fi = open(input_path, 'rb')
    except:
        ferr = open(os.path.join('temp', 'GF_error_status=-1'), 'w')
        ferr.close()
        return
    try:
        fo1 = open(output_path1,'wb', buffering=bufsize)
        fo2 = open(output_path2, 'wb', buffering=bufsize)
    except:
        ferr = open(os.path.join('temp', 'GF_error_status=-2'), 'w')
        ferr.close()
        return

    # variate initialization
    last_progress = 0
    headtext = fi.readline().decode('utf-8')
    progress = len(headtext)
    headtext = headtext.strip()
    headlist = headtext.split('\t')
    headlist = headlist[1:]
    strline = fi.readline().decode('utf-8')

    # write head text
    fo1.write((headtext + '\tP\n').encode('utf-8'))
    fo2.write((headtext + '\tP\tCorr\n').encode('utf-8'))

    # start time
    last_time = time.time()

    # main loop
    while strline != '':
        st = strline.find('\t') + 1
        si = strline[st:-1]
        si = si.split('\t')
        kmer_fre = []
        group1 = []
        group2 = []
        # logical filtering
        for i in range(len(si)):
            if not(headlist[i] in TI_dic):
                continue
            kmer_fre.append(float(si[i]))
            if si[i] == '0':
                group1.append(float(TI_dic[headlist[i]]))
            else:
                group2.append(float(TI_dic[headlist[i]]))
        (Zvalue, Pvalue) = stats.ranksums(group1, group2)
        if Pvalue < wicxon_p:
            fo1.write(strline[:-1].encode('utf-8'))
            fo1.write(('\t' + str(Pvalue) + '\n').encode('utf-8'))
        else:
            df = pd.DataFrame({'Kmer':kmer_fre, 'Trait':group1+group2})
            corr = df.corr('spearman')['Kmer'][1]
            if corr >= corr_value:
                fo2.write(strline[:-1].encode('utf-8'))
                fo2.write(('\t' + str(Pvalue) + '\t' + str(corr) + '\n').encode('utf-8'))

        progress += len(strline)     # update the progress
        if time.time() - last_time >= 1:    # output the progress file per second
            os.rename(os.path.join('temp', 'GF_progress' + str(Nprocess) + ' ' + str(last_progress)),
                      os.path.join('temp', 'GF_progress' + str(Nprocess) + ' ' + str(progress)))
            last_progress = progress
            last_time = time.time()
        strline = fi.readline().decode('utf-8')

    # progress 100%
    os.rename(os.path.join('temp', 'GF_progress' + str(Nprocess) + ' ' + str(last_progress)),
              os.path.join('temp', 'GF_progress' + str(Nprocess) + ' ' + str(progress)))

    # close file
    fi.close()
    fo1.close()
    fo2.close()

    # create ok-flag file
    fok = open(os.path.join('temp', str(Nprocess) + '_ok'), 'w')
    fok.close()

class GF_Thread(threading.Thread):
    def __init__(self, GF_Param):
        threading.Thread.__init__(self)
        self.status = 1
        self.loginfo = ''
        self.gm_result_path = GF_Param[0]
        self.gf_result_path = GF_Param[1]
        self.ass_l_value = GF_Param[2]
        self.p_value = GF_Param[3]
        self.ass_n_value = GF_Param[4]
        self.Number_of_A_Group = GF_Param[5]
        self.Number_of_B_Group = GF_Param[6]
        self.GroupA_Name = GF_Param[7]
        self.GroupB_Name = GF_Param[8]
        self.TI_dic = GF_Param[9]
        self.corr_value = GF_Param[10]
        self.Categorical_Mode = GF_Param[11]
        self.filesize = []
        self.file_path = []
        self.files_number = 0

    def run(self):
        filespath = os.listdir(self.gm_result_path)
        filespath.sort()
        for i in filespath:
            i = os.path.join(self.gm_result_path, i)
            if os.path.isdir(i) == False:
                self.filesize.append(os.path.getsize(i))
                self.file_path.append(i)
                self.files_number += 1

        if self.files_number > 256:
            self.status = -1
            self.loginfo = 'Input files number > 256.'
            return

        # create the temp folder
        if os.path.exists('temp') == False:
            os.mkdir('temp')

        for i in range(self.files_number):
            # the progress file
            fprogress = open(os.path.join('temp', 'GF_progress' + str(i) + ' 0'), 'w')
            fprogress.close()

        # the multiprocess runs
        categorical_param = [self.Number_of_A_Group, self.Number_of_B_Group, self.GroupA_Name, self.GroupB_Name,
                          self.TI_dic, self.ass_l_value, self.p_value, self.ass_n_value, 4096]
        continuous_param = [self.TI_dic, self.p_value, self.corr_value, 4096]
        if self.Categorical_Mode == True:
            self.jobs = [Process(target=Categorical_feature_filtering,
                            args=(i, self.file_path[i],
                                  os.path.join(self.gf_result_path, 'categorical_l_' + str(i) + '.txt'),
                                  os.path.join(self.gf_result_path, 'categorical_n_' + str(i) + '.txt'),
                                  categorical_param))
                    for i in range(self.files_number)]
        else:
            self.jobs = [Process(target=Continuous_feature_filtering,
                            args=(i, self.file_path[i],
                                  os.path.join(self.gf_result_path, 'continuous_l_' + str(i) + '.txt'),
                                  os.path.join(self.gf_result_path, 'continuous_n_' + str(i) + '.txt'),
                                  continuous_param))
                    for i in range(self.files_number)]
        self.status = 2
        _thread.start_new_thread(self.detective_process_ok, ())
        for j in self.jobs:
            j.start()
        for j in self.jobs:
            j.join()
        
        if self.detective_error() == 0:
            self.status = 0

        time.sleep(1)   # wait detective_process_ok to exit

    def detective_process_ok(self):
        while True:
            okflag = True
            for i in range(self.files_number):
                if not os.path.exists(os.path.join('temp', str(i) + '_ok')):
                    okflag = False
                    break
            if okflag == True or self.status != 2:
                try:
                    for j in self.jobs:
                        j.terminate()
                except:
                    pass
                return

    def detective_error(self):
        if os.path.exists(os.path.join('temp', 'GF_error_status=-1')):
            return -1
        elif os.path.exists(os.path.join('temp', 'GF_error_status=-2')):
            return -2
        return 0

    def detective_progress(self):
        if os.path.exists('temp') == True:
            progresslist = os.listdir('temp')
        else:
            progresslist = [0 for i in range(self.files_number - 1)]
        progress_now = [0 for i in range(self.files_number)]
        if self.status == 2:
            if len(progresslist) < self.files_number:
                self.status = -10
                self.loginfo = 'Missing progress files.'
                return progress_now
        else:
            return progress_now
        for i in progresslist:
            if i[:11] == 'GF_progress':
                colon_pos = i.find(' ')
                progress_now[int(i[11:colon_pos])] = int(i[colon_pos + 1:])
        return progress_now
