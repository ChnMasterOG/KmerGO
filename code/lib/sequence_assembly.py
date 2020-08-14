# coding = utf-8
# author: QiChen
# version: v2.2
# modification date: 2020/5/8
# assembly tool: cap3

import os, sys
import threading
import platform

class KA_Thread(threading.Thread):
    def __init__(self, KA_Param):
        threading.Thread.__init__(self)
        self.Parameters = KA_Param
        self.status = 1
        self.loginfo = ''
        self.cmdOption = ' -i 30  -j 31  -o 18  -s 300'

    def run(self):
        try:
            self.status = self.KA_GO()
        except:
            self.status = -999
        if self.status == -1:
            self.loginfo = 'The system is not supported.'
        elif self.status == -2:
            self.loginfo = 'Missing some feature files.'
        elif self.status == -3:
            self.loginfo = 'Missing the result folder.'
        elif self.status == -5:
            self.loginfo = 'Files format error.'
        elif self.status == -10:
            self.loginfo = 'No group-specific k-mer.'
        elif self.status == -999:
            self.loginfo = 'There are some errors when CAP3 is running.'

    def KA_GO(self):
        self.loginfo = 'CAP3 is working...'

        gf_path = self.Parameters[0]
        result_path = self.Parameters[1]
        system_platform = platform.system()

        # convert feature files to fasta
        try:
            foA = open(os.path.join(result_path, 'Afeatures.fa'), 'w', buffering=4096)
            foB = open(os.path.join(result_path, 'Bfeatures.fa'), 'w', buffering=4096)
        except:
            return -3
        NameA = ''
        NameB = ''
        countA = 0
        countB = 0
        filtering_mode = 0
        flist = []
        counter = 0
        for dir, folder, file in os.walk(gf_path):
            if dir == gf_path:
                flist = file
        for fpath in flist:
            counter += 1
            try:
                fi = open(os.path.join(gf_path, fpath), 'r')
            except:
                return -2
            if fpath[:11] == 'categorical':  # Categorical data
                if counter != 1 and filtering_mode == 1:
                    return -5
                filtering_mode = 0
            else:  # Continuous data
                if counter != 1 and filtering_mode == 0:
                    return -5
                filtering_mode = 1
            fi.readline()   # skip headtext
            s = fi.readline().strip()
            while s != '':
                if filtering_mode == 0:
                    label = s.split('\t')[-1]
                else:
                    label = 'Result'
                if NameA == '':
                    NameA = label
                elif NameB == '' and NameA != label:
                    NameB = label
                if label == NameA:
                    countA += 1
                    foA.write('>' + str(countA) + '\n')
                    foA.write(s[:s.find('\t')] + '\n')
                else:
                    countB += 1
                    foB.write('>' + str(countB) + '\n')
                    foB.write(s[:s.find('\t')] + '\n')
                s = fi.readline().strip()
            fi.close()
        foA.close()
        foB.close()

        if NameA != '':
            os.rename(os.path.join(result_path, 'Afeatures.fa'), os.path.join(result_path, NameA + '_specific.fa'))
        else:
            os.remove(os.path.join(result_path, 'Afeatures.fa'))
        if NameB != '':
            os.rename(os.path.join(result_path, 'Bfeatures.fa'), os.path.join(result_path, NameB + '_specific.fa'))
        else:
            os.remove(os.path.join(result_path, 'Bfeatures.fa'))
        if NameA == '' and NameB == '':
            return -10

        # using the cap3 tool
        if system_platform == 'Windows':
            cap3_command = r'.\bin\cap3.exe '
        elif system_platform == 'Linux':
            file_path = os.path.dirname(os.path.realpath(sys.argv[0]))
            cap3_command = file_path + '/bin/./cap3 '
        else:
            return -1  # the system is not supported

        if NameA != '':
            r1 = os.system(cap3_command + os.path.join(result_path, NameA + '_specific.fa') + self.cmdOption)
            os.rename(os.path.join(result_path, NameA + '_specific.fa'), os.path.join(result_path, NameA + '_specific_kmer.fa'))
        else:
            r1 = -1

        if NameB != '':
            r2 = os.system(cap3_command + os.path.join(result_path, NameB + '_specific.fa') + self.cmdOption)
            os.rename(os.path.join(result_path, NameB + '_specific.fa'), os.path.join(result_path, NameB + '_specific_kmer.fa'))
        else:
            r2 = -1
        
        if r1 != 0 and r2 != 0:
            return 0 / 0

        return 0

