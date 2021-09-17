# coding = utf-8
# author: QiChen
# version: v5.6
# modification date: 2020/10/5

import os, shutil
import time
import threading
import _thread
from multiprocessing import Process
from lib import loser_tree

Sparse_filter_threshold = 0.8

def distribute_kmer_nonuniformly(step, num, maxnum, out):
    if step == 4:
        return out
    if num < 4 * maxnum / 10:
        return distribute_kmer_nonuniformly(step + 1, num, 4 * maxnum / 10, out)
    elif num < 7 * maxnum / 10:
        out += 1 * (4 ** (3-step))
        return distribute_kmer_nonuniformly(step + 1, num - 4 * maxnum / 10, 3 * maxnum / 10, out)
    elif num < 9 * maxnum / 10:
        out += 2 * (4 ** (3-step))
        return distribute_kmer_nonuniformly(step + 1, num - 7 * maxnum / 10, 2 * maxnum / 10, out)
    else:
        out += 3 * (4 ** (3-step))
        return distribute_kmer_nonuniformly(step + 1, num - 9 * maxnum / 10, maxnum / 10, out)

def get_Son_Matrix(Nprocess, NEXTprocess, param):
    gm_result_path = param[0]
    beacon_path_list = param[1]
    path_list = param[2]
    Klen = param[3]
    Number_of_Group = param[4]
    head_list = param[5]
    A_Name = param[6]
    A_Number = param[7]
    TI_dic = param[8]
    KofZ = param[9]
    fre_sum = param[10]
    bufsize = param[11]

    # Memory exchange efficiency
    if A_Name is not None:
        no_zero_counter_thrA = (1 - Sparse_filter_threshold) * A_Number
        if A_Name == '':
            no_zero_counter_thrB = 1
        else:
            no_zero_counter_thrB = (1 - Sparse_filter_threshold) * (Number_of_Group - A_Number)
    else:
        no_zero_counter_thrA = 1
        no_zero_counter_thrB = 1
    next_Nprocess = Nprocess + NEXTprocess
    Klen_and_1 = Klen + 1
    KofZ_t0 = KofZ + '\t0'
    zero_matrix = ['0' for i in range(Number_of_Group)]

    ##########################################
    # STEP1: locate to the interruputed_kmer #
    ##########################################
    beacon_point = [0 for i in range(len(beacon_path_list))]
    beacon_sum = [0 for i in range(len(beacon_path_list))]
    for i in range(len(beacon_path_list)):
        try:
            f_beacon = open(beacon_path_list[i], 'r')
        except:
            ferr = open(os.path.join('temp', 'GM_error_status=-1'), 'w')
            ferr.close()
            return
        text = f_beacon.readline().strip()
        while text != 'beacon:':
            text = f_beacon.readline().strip()
        for j in range(Nprocess + 1):
            beacon_point[i] = int(f_beacon.readline().strip())
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
            ferr = open(os.path.join('temp', 'GM_error_status=-1'), 'w')
            ferr.close()
            return
        try:
            f[i].seek(beacon_point[i], 0)
        except:
            ferr = open(os.path.join('temp', 'GM_error_status=-3'), 'w')
            ferr.close()
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

    # open the output file
    try:
        fout = open(os.path.join(gm_result_path, 'son_matrix_' + str(Nprocess) + '.txt'), 'wb', buffering=bufsize)
    except:
        ferr = open(os.path.join('temp', 'GM_error_status=-2'), 'w')
        ferrc.close()
        return

    # write head text
    fout.write(('\t'.join(head_list)).encode('utf-8'))
    head_list = head_list[1:]

    # read each lines from files
    s = [f[i].readline().decode('utf-8') for i in range(Number_of_Group)]

    # preliminary judgment
    for i in range(len(s)):
        if s[i] == '' or s[i][:4] >= BEACON_prefix[next_Nprocess]:
            s[i] = KofZ_t0
            END_flag += 1
        else:
            progress += len(s[i])    # update progress
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
            temp_s_min_index = s[min_index]
            temp_min_kmer = temp_s_min_index[:Klen]
            if min_kmer == temp_min_kmer:  # new k-mer is equal to last k-mer
                wline[min_index + 1] = ('%.4f' % (int(temp_s_min_index[Klen_and_1:]) / fre_sum[min_index]))
                if A_Name is not None:
                    if A_Name == '':
                        no_zero_counter1 += 1
                    else:
                        if TI_dic[head_list[min_index]] == A_Name:
                            no_zero_counter1 += 1
                        else:
                            no_zero_counter2 += 1
                else:
                    no_zero_counter1 += 1
            else:
                if no_zero_counter1 >= no_zero_counter_thrA or no_zero_counter2 >= no_zero_counter_thrB:
                    fout.write(('\n' + '\t'.join(wline)).encode('utf-8'))  # write last wline
                no_zero_counter1 = no_zero_counter2 = 0
                wline = [temp_min_kmer]
                wline.extend(zero_matrix)
                wline[min_index + 1] = ('%.4f' % (int(temp_s_min_index[Klen_and_1:]) / fre_sum[min_index]))
                if A_Name is not None:
                    if A_Name == '':
                        no_zero_counter1 += 1
                    else:
                        if TI_dic[head_list[min_index]] == A_Name:
                            no_zero_counter1 += 1
                        else:
                            no_zero_counter2 += 1
                else:
                    no_zero_counter1 += 1
            min_kmer = temp_min_kmer
            s[min_index] = f[min_index].readline().decode('utf-8')
            temp_s_min_index = s[min_index]
            if time.time() - last_time >= 5:    # output the progress file every 5 seconds
                os.rename(os.path.join('temp', 'GM_progress' + str(Nprocess) + ' ' + str(last_progress)),
                          os.path.join('temp', 'GM_progress' + str(Nprocess) + ' ' + str(progress)))
                last_progress = progress
                last_time = time.time()
            if temp_s_min_index == '' or temp_s_min_index[:4] >= BEACON_prefix[next_Nprocess]:
                s[min_index] = KofZ_t0
                END_flag += 1
            else:
                progress += len(temp_s_min_index)    # update progress
            if END_flag == Number_of_Group:
                if no_zero_counter1 >= no_zero_counter_thrA or no_zero_counter2 >= no_zero_counter_thrB:
                    fout.write(('\n' + '\t'.join(wline)).encode('utf-8'))  # write last wline
                fout.write(b'\n')  # write END line
                break
            losertree.set_kmerlist(min_index, s[min_index])
            losertree.adjust(min_index)
    except:
        ferr = open(os.path.join('temp', 'GM_error_status=-10'), 'w')
        ferr.close()
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

    # create ok-flag file
    fok = open(os.path.join('temp', str(Nprocess) + '_ok'), 'w')
    fok.close()

class GM_Thread(threading.Thread):
    def __init__(self, GM_Param):
        threading.Thread.__init__(self)
        self.status = 1
        self.loginfo = ''
        self.filesize = []
        self.path_list = []
        self.beacon_path_list = []
        self.kmc_result_path = GM_Param[0]
        self.gm_result_path = GM_Param[1]
        self.process_number = GM_Param[2]
        self.A_Name = GM_Param[3]
        self.A_Number = GM_Param[4]
        self.TI_dic = GM_Param[5]
        self.KofZ = ''
        self.Klen = -1

    def run(self):
        self.path_list = []
        self.beacon_path_list = []
        self.beacon_block_list = []
        self.head_list = ['k-mer']
        flist = []
        try:
            for dir, folder, file in os.walk(self.kmc_result_path):
                if dir == self.kmc_result_path:
                    flist = file
            for i in list(flist):
                if i[-11:] == '_beacon.txt':
                    self.beacon_path_list.append(os.path.join(self.kmc_result_path, i))
                    if i[:-11] + '.txt' not in flist:
                        self.status = -70
                        self.loginfo = 'Missing some k-mer counting files.'
                        return
                    self.head_list.append(i[:-11])
                    self.path_list.append(os.path.join(self.kmc_result_path, i[:-11] + '.txt'))
                    self.filesize.append(os.path.getsize(os.path.join(self.kmc_result_path, i[:-11] + '.txt')))
                    # read kmc files to set Klen
                    ftmp = open(os.path.join(self.kmc_result_path, i[:-11] + '.txt'), 'r')
                    firstline = ftmp.readline().strip()
                    firstline = firstline.split('\t')
                    if len(firstline[0]) != 0:
                        if self.Klen != -1 and self.Klen != len(firstline[0]):
                            ftmp.close()
                            self.status = -80
                            self.loginfo = 'K value error.'
                            return
                        else:
                            self.Klen = len(firstline[0])
                    else:   # fill beacon
                        f_beacon_tmp = open(os.path.join(self.kmc_result_path, i), 'w')
                        f_beacon_tmp.write('beacon:\n')
                        for k in range(256):
                            f_beacon_tmp.write('0\n')
                        f_beacon_tmp.write('sum:\n0\n')
                        f_beacon_tmp.close()
                    ftmp.close()
        except:
            self.status = -3
            self.loginfo = 'Can not open some KMC result files.'
            return

        self.KofZ = 'Z' * self.Klen

        if len(self.beacon_path_list) != len(self.path_list):
            self.status = -15
            self.loginfo = 'Files\' number is wrong.'
            return

        # create the temp folder
        if os.path.exists('temp') == False:
            os.mkdir('temp')
        for i in range(self.process_number):
            if self.process_number >= 40:
                self.beacon_block_list.append(round(i * 256 / self.process_number))
            else:
                self.beacon_block_list.append(distribute_kmer_nonuniformly(0, round(i * 10000 / self.process_number), 10000, 0))
            # the progress file
            fprogress = open(os.path.join('temp', 'GM_progress' + str(self.beacon_block_list[i]) + ' 0'), 'w')
            fprogress.close()
        self.beacon_block_list.append(256)

        # get each block size
        fre_sum = []
        beacon_point = [[0 for j in range(self.process_number + 1)] for i in range(len(self.path_list))]
        block_size = [[0 for j in range(len(self.path_list))] for i in range(self.process_number)]
        for i in range(len(self.path_list)):
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
                beacon_point[i][j] = int(text)
            try:
                beacon_point[i][self.process_number] = int(self.filesize[i])
            except:
                self.status = -99
                self.loginfo = 'K vaule or Cs value is wrong.'
                return
            for j in range(self.process_number):
                block_size[j][i] = beacon_point[i][j + 1] - beacon_point[i][j]
            f_beacon.seek(0, 0)
            text = f_beacon.readline().strip()
            while text != 'sum:':
                text = f_beacon.readline().strip()
            fre_sum.append(int(f_beacon.readline().strip()))
            f_beacon.close()
        self.block_size = [sum(n) for n in block_size]

        # normalization coefficient
        tmp = -5
        for i in range(len(self.path_list)):
            tmp = max(tmp, len(str(fre_sum[i])) - 4)
        for i in range(len(self.path_list)):
            if tmp > 0:
                fre_sum[i] = fre_sum[i] / (10 ** tmp)
            elif tmp < 0:
                fre_sum[i] = fre_sum[i] * (10 ** -tmp)

        # output the normalization coefficients
        f_NC = open('normalization_coefficients.txt', 'w')
        for i in self.head_list:
            if i == 'k-mer':
                f_NC.write('sample')
            else:
                f_NC.write('\t' + i)
        f_NC.write('\nnormalization coefficient')
        for i in fre_sum:
            f_NC.write('\t' + str(i))
        f_NC.write('\n')
        f_NC.close()

        # the multiprocess runs
        process_param = [self.gm_result_path, self.beacon_path_list, self.path_list, self.Klen,
                         len(self.path_list), self.head_list, self.A_Name, self.A_Number, self.TI_dic, self.KofZ,
                         fre_sum, 4096]
        self.jobs = [Process(target = get_Son_Matrix,
                             args = (self.beacon_block_list[i], self.beacon_block_list[i + 1] - self.beacon_block_list[i],
                                     process_param))
                for i in range(self.process_number)
        ]
        self.status = 2
        _thread.start_new_thread(self.detective_process_ok, ())
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

        time.sleep(1)   # wait detective_process_ok to exit
        
    def detective_process_ok(self):
        while True:
            okflag = True
            for i in range(self.process_number):
                if not os.path.exists(os.path.join('temp', str(self.beacon_block_list[i]) + '_ok')):
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
                self.status = -11
                self.loginfo = 'Missing progress files.'
                return progress_now
        else:
            return progress_now
        for i in progresslist:
            if i[:11] == 'GM_progress':
                colon_pos = i.find(' ')
                progress_now[self.beacon_block_list.index(int(i[11:colon_pos]))] = int(i[colon_pos + 1:])
        return progress_now
