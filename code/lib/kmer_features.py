# coding = utf-8
# author: QiChen
# version: v4.5
# modification date: 2019/12/5

import os
import time
import threading
from scipy import stats
import numpy as np
from multiprocessing import Process
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix

def feature_filtering(Nprocess, input_path, output_path1, output_path2, param):
    Number_of_A_Group = param[0]
    Number_of_B_Group = param[1]
    ass_l = param[2]
    wicxon_p = param[3]
    ass_n = param[4]
    bufsize = param[5]

    try:
        fi = open(input_path, 'r')
    except:
        open(os.path.join('temp', 'GF_error_status=-1'), 'w')
        return
    try:
        fo1 = open(output_path1,'w', buffering=bufsize)
        fo2 = open(output_path2, 'w', buffering=bufsize)
    except:
        open(os.path.join('temp', 'GF_error_status=-2'), 'w')
        return

    # variate initialization
    progress = 0
    last_progress = 0
    strline = fi.readline()

    # start time
    last_time = time.time()

    # main loop
    while strline != '':
        st = strline.find('\t') + 1
        si = strline[st:-1]
        si = si.split('\t')
        so = strline[:st]
        tp = 0
        fn = 0
        group1 = []
        group2 = []
        label_Y = []
        # logical filtering
        for i in range(len(si)):
            nor_value = float(si[i])
            so += str(nor_value) + '\t'
            if i < Number_of_A_Group:
                group1.append(nor_value)
                label_Y.append(0)
                if si[i] != '0':
                    tp += 1
            else:
                group2.append(nor_value)
                label_Y.append(1)
                if si[i] != '0':
                    fn += 1
        so = so[:-1]
        ASStrue = tp / Number_of_A_Group + (Number_of_B_Group - fn) / Number_of_B_Group
        ASStrue /= 2
        if ASStrue > 1 - ASStrue:
            label = 'A'
        else:
            ASStrue = 1 - ASStrue
            label = 'B'
        if ASStrue >= ass_l:
            fo1.write(so)
            fo1.write('\tLogical ASS:' + str(ASStrue) + '\tlabel:' + str(label) + '\n')
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
                        label = 'A'
                    else:
                        label = 'B'
                    fo2.write(so)
                    fo2.write('\tLogical_ASS:' + str(ASStrue) + '\tWilcoxon_p:' + str(Pvalue) + '\tNumerical_ASS:' +
                              str(cc) + '\tlabel:' + str(label) + '\n')
        progress += len(strline)     # update the progress
        if time.time() - last_time >= 1:    # output the progress file per second
            os.rename(os.path.join('temp', 'GF_progress' + str(Nprocess) + ' ' + str(last_progress)),
                      os.path.join('temp', 'GF_progress' + str(Nprocess) + ' ' + str(progress)))
            last_progress = progress
            last_time = time.time()
        strline = fi.readline()

    # progress 100%
    os.rename(os.path.join('temp', 'GF_progress' + str(Nprocess) + ' ' + str(last_progress)),
              os.path.join('temp', 'GF_progress' + str(Nprocess) + ' ' + str(progress)))

    # close file
    fi.close()
    fo1.close()
    fo2.close()

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
            open(os.path.join('temp', 'GF_progress' + str(i) + ' 0'), 'w')

        # the multiprocess runs
        process_param = [self.Number_of_A_Group, self.Number_of_B_Group, self.ass_l_value,
                         self.p_value, self.ass_n_value, 4096]
        self.jobs = [Process(target=feature_filtering,
                        args=(i, self.file_path[i],
                              os.path.join(self.gf_result_path, 'logical_' + str(i) + '.txt'),
                              os.path.join(self.gf_result_path, 'numeric_' + str(i) + '.txt'),
                              process_param))
                for i in range(self.files_number)]
        self.status = 2
        for j in self.jobs:
            j.start()
        for j in self.jobs:
            j.join()
        
        if self.detective_error() == 0:
            self.status = 0

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
