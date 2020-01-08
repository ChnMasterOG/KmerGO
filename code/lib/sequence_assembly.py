# coding = utf-8
# author: QiChen
# version: v2.0
# modification date: 2019/12/6
# assembly tool: cap3

import os
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
        countA = 0
        countB = 0
        block_number = int(len(os.listdir(gf_path)) / 2)
        for i in range(block_number):
            try:
                fi = open(os.path.join(gf_path, 'logical_' + str(i) + '.txt'), 'r')
            except:
                return -2
            s = fi.readline()
            while s != '':
                if s[-2] == 'A':
                    countA += 1
                    foA.write('>' + str(countA) + '\n')
                    foA.write(s[:s.find('\t')] + '\n')
                else:
                    countB += 1
                    foB.write('>' + str(countB) + '\n')
                    foB.write(s[:s.find('\t')] + '\n')
                s = fi.readline()
            fi.close()
            try:
                fi = open(os.path.join(gf_path, 'numeric_' + str(i) + '.txt'), 'r')
            except:
                return -2
            s = fi.readline()
            while s != '':
                if s[-2] == 'H':
                    countA += 1
                    foA.write('>' + str(countA) + '\n')
                    foA.write(s[:s.find('\t')] + '\n')
                else:
                    countB += 1
                    foB.write('>' + str(countB) + '\n')
                    foB.write(s[:s.find('\t')] + '\n')
                s = fi.readline()
            fi.close()
        foA.close()
        foB.close()

        # using the cap3 tool
        if system_platform == 'Windows':
            cap3_command = r'.\bin\cap3.exe '
        elif system_platform == 'Linux':
            cap3_command = './bin/cap3 '
        else:
            return -1  # the system is not supported

        r1 = os.system(cap3_command + os.path.join(result_path, 'Afeatures.fa') + self.cmdOption)

        r2 = os.system(cap3_command + os.path.join(result_path, 'Bfeatures.fa') + self.cmdOption)
        
        if r1 != 0 and r2 != 0:
            return 0 / 0

        return 0

