# coding = utf-8
# author: QiChen
# version: v2.6
# modification date: 2020/1/14

import os
import subprocess
import threading
import platform

FASTA_suffix = ['.fa', '.fna', '.fasta', '.fa.gz', '.fna.gz', '.fasta.gz']
FASTQ_suffix = ['.fq', '.fastq', '.fq.gz', '.fastq.gz']

class KMC_Thread(threading.Thread):
    def __init__(self, KMC_Param):
        threading.Thread.__init__(self)
        self.Parameters = KMC_Param
        self.loginfo = ''
        self.status = 1

    def run(self):
        try:
            self.status = self.KMC_GO()
        except:
            self.status = -999
        if self.status == -1:
            self.loginfo = 'The system is not supported.'
        elif self.status == -2:
            self.loginfo = 'There are no FASTA/Q files under the selected folder.'
        elif self.status == -999:
            self.loginfo = 'There are some errors when KMC is running.'

    def KMC_GO(self):
        self.loginfo = 'KMC working is ready...'

        system_platform = platform.system()
        k_value = self.Parameters[0]
        ci_value = self.Parameters[1]
        cs_value = self.Parameters[2]
        work_dir = [self.Parameters[3], self.Parameters[4]]
        output_dir = self.Parameters[5]
        key = ['A', 'B']
        flist = [[], []]
        typelist = []

        if system_platform == 'Windows':
            kmc_command = r'.\bin\kmc.exe'
            kmc_tools_command = r'.\bin\kmc_tools.exe'
            kmc_dump_command = r'.\bin\kmc_dump.exe'
        elif system_platform == 'Linux':
            kmc_command = './bin/kmc'
            kmc_tools_command = './bin/kmc_tools'
            kmc_dump_command = './bin/kmc_dump'
        else:
            return -1   # the system is not supported

        for j in range(2):  # A/B group
            for dir, folder, file in os.walk(work_dir[j]):
                if dir == work_dir[j]:
                    flist[j] = file
            flist[j].sort()
            for i in list(flist[j]):
                suffix = i[i.find('.'):]
                if suffix.lower() in FASTA_suffix:
                    typelist.append(' -fm ')
                    flist[j][flist[j].index(i)] = os.path.join(work_dir[j], i)
                elif suffix.lower() in FASTQ_suffix:
                    typelist.append(' -fq ')
                    flist[j][flist[j].index(i)] = os.path.join(work_dir[j], i)
                else:
                    flist[j].remove(i)
            if len(flist[j]) == 0:
                return -2  # no files in this directory

        for j in range(2):  # A/B group
            countings = 0
            self.loginfo = key[j] + ' group 0/' + str(len(flist[j]))
            for i in range(len(flist[j])):
                countings += 1
                cmd = kmc_command + ' -k' + str(k_value) + ' -ci' + str(ci_value) + ' -cs' + str(cs_value)\
                      + typelist[i] + flist[j][i] + ' ' + os.path.join(work_dir[j], key[j] + str(countings))\
                      + ' ' + work_dir[j]
                r = subprocess.call(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
                if r != 0:
                    return 0 / 0
                cmd = kmc_tools_command + ' transform ' + os.path.join(work_dir[j], key[j] + str(countings))\
                      + ' sort ' + os.path.join(work_dir[j], key[j] + str(countings) + '_sort')
                r = subprocess.call(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
                if r != 0:
                    return 0 / 0
                cmd = kmc_dump_command + ' -cs' + str(cs_value) + ' '\
                      + os.path.join(work_dir[j], key[j] + str(countings) + '_sort') + ' '\
                      + os.path.join(output_dir, key[j] + str(countings) + '.txt') + ' '\
                      + os.path.join(output_dir, key[j]+ str(countings) + '_beacon.txt')
                r = subprocess.call(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
                if r != 0:
                    return 0 / 0
                os.remove(os.path.join(work_dir[j], key[j] + str(countings) + '.kmc_pre'))
                os.remove(os.path.join(work_dir[j], key[j] + str(countings) + '.kmc_suf'))
                os.remove(os.path.join(work_dir[j], key[j] + str(countings) + '_sort.kmc_pre'))
                os.remove(os.path.join(work_dir[j], key[j] + str(countings) + '_sort.kmc_suf'))
                self.loginfo = key[j] + ' group ' + str(countings) + '/' +  str(len(flist[j]))

        return 0
