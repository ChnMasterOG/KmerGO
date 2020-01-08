# coding = utf-8
# author: QiChen
# version: v5.1
# modification date: 2019/12/12

import os, shutil
import time
import threading
from multiprocessing import Process
from lib import loser_tree

Sparse_filter_threshold = 0.8

def get_Son_Matrix(Nprocess, NEXTprocess, param):
    gm_result_path = param[0]
    beacon_path_list = param[1]
    path_list = param[2]
    Klen = param[3]
    fre_max_len = param[4]
    Number_of_A_Group = param[5]
    Number_of_B_Group = param[6]
    KofZ = param[7]
    fre_sum = param[8]
    bufsize = param[9]

    # Memory exchange efficiency
    Number_of_Group = Number_of_A_Group + Number_of_B_Group
    no_zero_counter_thrA = (1 - Sparse_filter_threshold) * Number_of_A_Group
    no_zero_counter_thrB = (1 - Sparse_filter_threshold) * Number_of_B_Group
    next_Nprocess = Nprocess + NEXTprocess
    line_length = Klen + fre_max_len + 2
    Klen_and_1 = Klen + 1
    Klen_and_2 = Klen_and_1 + 1
    KofZ_t0 = KofZ + '\t0'
    zero_matrix = ['0' for i in range(Number_of_Group)]

    ##########################################
    # STEP1: locate to the interruputed_kmer #
    ##########################################
    beacon_line = [0 for i in range(len(beacon_path_list))]
    beacon_sum = [0 for i in range(len(beacon_path_list))]
    for i in range(len(beacon_path_list)):
        try:
            f_beacon = open(beacon_path_list[i], 'r')
        except:
            open(os.path.join('temp', 'GM_error_status=-1'), 'w')
            return
        text = f_beacon.readline().strip()
        while text != 'beacon:':
            text = f_beacon.readline().strip()
        for j in range(Nprocess + 1):
            beacon_line[i] = int(f_beacon.readline().strip())
        while text != 'sum:':
            text = f_beacon.readline().strip()
        beacon_sum[i] = int(f_beacon.readline().strip())
        f_beacon.close()
    f = []
    for i in range(len(path_list)):
        try:
            f.append(open(path_list[i], 'rb'))
        except:
            for j in range(i):
                f[i].close()
            open(os.path.join('temp', 'GM_error_status=-1'), 'w')
            return
        try:
            f[i].seek(line_length * beacon_line[i], 0)
        except:
            open(os.path.join('temp', 'GM_error_status=-3'), 'w')
            return

    #####################################
    # STEP2: set block to build a matrix#
    #####################################
    # constant variate
    BEACON_prefix = ['' for i in range(256 + 1)]
    for i in range(256):
        temp = 256
        for j in range(4):
            temp /= 4
            if (i // temp) % 4 == 0:
                BEACON_prefix[i] += 'A'
            elif (i // temp) % 4 == 1:
                BEACON_prefix[i] += 'C'
            elif (i // temp) % 4 == 2:
                BEACON_prefix[i] += 'G'
            else:
                BEACON_prefix[i] += 'T'
    for i in range(4):
        BEACON_prefix[256] += 'Z'

    # variate initialization
    END_flag = 0
    no_zero_counter1 = 0
    no_zero_counter2 = 0
    min_kmer = ''
    wline = ['']
    progress = 0
    last_progress = 0
    const_strlen_list = [Klen_and_1 + 2 * i for i in range(Number_of_Group)]
    wline_strlen_list = const_strlen_list

    # open the output file
    try:
        fout = open(os.path.join(gm_result_path, 'son_matrix_' + str(Nprocess) + '.txt'), 'wb', buffering=bufsize)
    except:
        open(os.path.join('temp', 'GM_error_status=-2'), 'w')
        return

    # read each lines from files
    s = [f[i].read(line_length).decode('utf-8') for i in range(Number_of_Group)]

    # preliminary judgment
    for i in range(len(s)):
        if s[i] == '' or s[i][:4] >= BEACON_prefix[next_Nprocess]:
            s[i] = KofZ_t0
            END_flag += 1
    if END_flag == Number_of_Group:
        # close files
        for i in range(len(path_list)):
            f[i].close()
        fout.close()
        return

    # build a loser tree
    losertree = loser_tree.LoserTree(s, Klen, Number_of_Group)
    losertree.setup_param()
    losertree.build()

    # start time
    last_time = time.time()
    
    try:
        # main loop
        while True:
            min_index = losertree.getmin()
            min_index_and_1 = min_index + 1
            if min_kmer == s[min_index][:Klen]:  # new k-mer is equal to last k-mer
                wline[min_index_and_1] = fre_str = ('%.4f' % (int(s[min_index][Klen_and_1:]) / fre_sum[min_index]))
                if min_index < Number_of_A_Group:
                    no_zero_counter1 += 1
                else:
                    no_zero_counter2 += 1
            else:
                if no_zero_counter1 >= no_zero_counter_thrA or no_zero_counter2 >= no_zero_counter_thrB:
                    fout.write(('\t'.join(wline) + '\n').encode('utf-8'))  # write last wline
                no_zero_counter1 = no_zero_counter2 = 0
                wline_strlen_list = const_strlen_list
                wline = [s[min_index][:Klen]]
                wline.extend(zero_matrix)
                wline[min_index_and_1] = ('%.4f' % (int(s[min_index][Klen_and_1:]) / fre_sum[min_index]))
                if min_index < Number_of_A_Group:
                    no_zero_counter1 += 1
                else:
                    no_zero_counter2 += 1
            min_kmer = s[min_index][:Klen]
            s[min_index] = f[min_index].read(line_length).decode('utf-8')
            progress += line_length    # update progress
            if time.time() - last_time >= 5:    # output the progress file every 5 seconds
                os.rename(os.path.join('temp', 'GM_progress' + str(Nprocess) + ' ' + str(last_progress)),
                          os.path.join('temp', 'GM_progress' + str(Nprocess) + ' ' + str(progress)))
                last_progress = progress
                last_time = time.time()
            if s[min_index] == '' or s[min_index][:4] >= BEACON_prefix[next_Nprocess]:
                s[min_index] = KofZ_t0
                END_flag += 1
            if END_flag == Number_of_Group:
                if no_zero_counter1 >= no_zero_counter_thrA or no_zero_counter2 >= no_zero_counter_thrB:
                    fout.write(('\t'.join(wline) + '\n').encode('utf-8'))  # write last wline
                break
            losertree.set_kmerlist(min_index, s[min_index])
            losertree.adjust(min_index)
    except:
        open(os.path.join('temp', 'GM_error_status=-10'), 'w')
        return
    
    # free the memory
    losertree.free_mem()

    # progress 100%
    os.rename(os.path.join('temp', 'GM_progress' + str(Nprocess) + ' ' + str(last_progress)),
              os.path.join('temp', 'GM_progress' + str(Nprocess) + ' ' + str(progress)))

    # close files
    for i in range(len(path_list)):
        f[i].close()
    fout.close()

class GM_Thread(threading.Thread):
    def __init__(self, GM_Param):
        threading.Thread.__init__(self)
        self.status = 1
        self.loginfo = ''
        self.filesize = []
        self.Number_of_A_Group = 0
        self.Number_of_B_Group = 0
        self.path_list = []
        self.beacon_path_list = []
        self.kmc_result_path = GM_Param[0]
        self.gm_result_path = GM_Param[1]
        self.process_number = GM_Param[2]
        self.Klen = GM_Param[3]
        cs = GM_Param[4]
        self.Number_of_A_Group = GM_Param[5]
        self.Number_of_B_Group = GM_Param[6]
        cs_len = 0
        while cs > 0:
            cs_len += 1
            cs = int(cs / 10)
        self.fre_max_len = cs_len
        self.KofZ = ''
        for i in range(self.Klen):
            self.KofZ += 'Z'
        self.Number_of_Group = self.Number_of_A_Group + self.Number_of_B_Group

    def run(self):
        self.path_list = []
        self.beacon_path_list = []
        self.beacon_block_list = []
        try:
            for i in range(1, self.Number_of_A_Group + 1):
                self.path_list.append(os.path.join(self.kmc_result_path, 'A' + str(i) + '.txt'))
                self.beacon_path_list.append(os.path.join(self.kmc_result_path, 'A' + str(i) + '_beacon.txt'))
                self.filesize.append(os.path.getsize(os.path.join(self.kmc_result_path, 'A' + str(i) + '.txt')))
            for i in range(1, self.Number_of_B_Group + 1):
                self.path_list.append(os.path.join(self.kmc_result_path, 'B' + str(i) + '.txt'))
                self.beacon_path_list.append(os.path.join(self.kmc_result_path, 'B' + str(i) + '_beacon.txt'))
                self.filesize.append(os.path.getsize(os.path.join(self.kmc_result_path, 'B' + str(i) + '.txt')))
        except:
            self.status = -3
            self.loginfo = 'Can not open some KMC result files.'
            return

        # create the temp folder
        if os.path.exists('temp') == False:
            os.mkdir('temp')
        for i in range(self.process_number):
            self.beacon_block_list.append(round(i * 256 / self.process_number))
            # the progress file
            open(os.path.join('temp', 'GM_progress' + str(self.beacon_block_list[i]) + ' 0'), 'w')
        self.beacon_block_list.append(256)

        # get each block size
        fre_sum = []
        beacon_line = [[0 for j in range(self.process_number + 1)] for i in range(self.Number_of_Group)]
        block_size = [[0 for j in range(self.Number_of_Group)] for i in range(self.process_number)]
        for i in range(self.Number_of_Group):
            try:
                f_beacon = open(self.beacon_path_list[i], 'r')
            except:
                self.status = -1
                self.loginfo = 'Can not open the K-mer beacon files.'
                return
            text = f_beacon.readline().strip()
            while text != 'beacon:':
                text = f_beacon.readline().strip()
            for j in range(self.process_number):
                if j != 0:
                    for k in range(self.beacon_block_list[j] - self.beacon_block_list[j - 1] - 1):
                        f_beacon.readline().strip()
                text = f_beacon.readline().strip()
                beacon_line[i][j] = int(text)
            try:
                beacon_line[i][self.process_number] = int(self.filesize[i] / (self.Klen + self.fre_max_len + 2))
            except:
                self.status = -99
                self.loginfo = 'K vaule or Cs value is wrong.'
                return
            for j in range(self.process_number):
                block_size[j][i] = (beacon_line[i][j + 1] - beacon_line[i][j]) * (self.Klen + self.fre_max_len + 2)
            f_beacon.seek(0, 0)
            text = f_beacon.readline().strip()
            while text != 'sum:':
                text = f_beacon.readline().strip()
            fre_sum.append(int(f_beacon.readline().strip()))
            f_beacon.close()
        self.block_size = [sum(n) for n in block_size]

        # normalization coefficient
        tmp = -5
        for i in range(self.Number_of_Group):
            tmp = max(tmp, len(str(fre_sum[i])) - 4)
        for i in range(self.Number_of_Group):
            if tmp > 0:
                fre_sum[i] = fre_sum[i] / (10 ** tmp)
            elif tmp < 0:
                fre_sum[i] = fre_sum[i] * (10 ** -tmp)

        # the multiprocess runs
        process_param = [self.gm_result_path, self.beacon_path_list, self.path_list, self.Klen, self.fre_max_len,
                         self.Number_of_A_Group, self.Number_of_B_Group, self.KofZ, fre_sum, 4096]
        self.jobs = [Process(target = get_Son_Matrix,
                             args = (self.beacon_block_list[i], self.beacon_block_list[i + 1] - self.beacon_block_list[i],
                                     process_param))
                for i in range(self.process_number)
        ]
        self.status = 2
        for j in self.jobs:
            j.start()
        for j in self.jobs:
            j.join()

        # rename son matrix
        j = 0
        self.beacon_block_list = self.beacon_block_list[:-1]
        for i in self.beacon_block_list:
            try:
                shutil.move(os.path.join(self.gm_result_path, 'son_matrix_' + str(i) + '.txt'),
                            os.path.join(self.gm_result_path, 'son_matrix_' + str(j) + '.txt'))
            except:
                pass
            j += 1

        if self.detective_error() == 0:
            self.status = 0

    def detective_error(self):
        if os.path.exists(os.path.join('temp', 'GM_error_status=-1')):
            return -1
        elif os.path.exists(os.path.join('temp', 'GM_error_status=-2')):
            return -2
        elif os.path.exists(os.path.join('temp', 'GM_error_status=-3')):
            return -3
        elif os.path.exists(os.path.join('temp', 'GM_error_status=-10')):
            return -10
        return 0

    def detective_progress(self):
        if os.path.exists('temp') == True:
            progresslist = os.listdir('temp')
        else:
            progresslist = [0 for i in range(self.process_number - 1)]
        progress_now = [0 for i in range(self.process_number)]
        if self.status == 2:
            if len(progresslist) < self.process_number:
                self.status = -10
                self.loginfo = 'Missing progress files.'
                return progress_now
        else:
            return progress_now
        for i in progresslist:
            if i[:11] == 'GM_progress':
                colon_pos = i.find(' ')
                progress_now[self.beacon_block_list.index(int(i[11:colon_pos]))] = int(i[colon_pos + 1:])
        return progress_now
