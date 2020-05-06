# coding = utf-8
# author: QiChen
# version: v3.0
# modification date: 2020/4/21

import os, sys
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
        work_dir = self.Parameters[3]
        output_dir = self.Parameters[4]
        flist = []
        typelist = []

        if system_platform == 'Windows':
            kmc_command = r'.\bin\kmc.exe'
            kmc_tools_command = r'.\bin\kmc_tools.exe'
            kmc_dump_command = r'.\bin\kmc_dump.exe'
        elif system_platform == 'Linux':
            file_path = os.path.dirname(os.path.realpath(sys.argv[0]))
            kmc_command = file_path + '/bin/./kmc'
            kmc_tools_command = file_path + '/bin/./kmc_tools'
            kmc_dump_command = file_path + '/bin/./kmc_dump'
        else:
            return -1   # the system is not supported

        for dir, folder, file in os.walk(work_dir):
            if dir == work_dir:
                flist = file
        flist.sort()
        for i in list(flist):
            suffix = i[i.find('.'):]
            if suffix.lower() in FASTA_suffix:
                typelist.append(' -fm ')
                flist[flist.index(i)] = [i[:i.find('.')], os.path.join(work_dir, i)]
            elif suffix.lower() in FASTQ_suffix:
                typelist.append(' -fq ')
                flist[flist.index(i)] = [i[:i.find('.')], os.path.join(work_dir, i)]
            else:
                flist.remove(i)
            if len(flist) == 0:
                return -2  # no files in this directory

        countings = 0
        self.loginfo = 'Number 0/' + str(len(flist))
        for i in range(len(flist)):
            countings += 1
            cmd = kmc_command + ' -k' + str(k_value) + ' -ci' + str(ci_value) + ' -cs' + str(cs_value)\
                  + typelist[i] + flist[i][1] + ' ' + os.path.join(work_dir, flist[i][0])\
                  + ' ' + work_dir
            r = subprocess.call(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
            if r != 0:
                return 0 / 0
            cmd = kmc_tools_command + ' transform ' + os.path.join(work_dir, flist[i][0])\
                  + ' sort ' + os.path.join(work_dir, flist[i][0] + '_sort')
            r = subprocess.call(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
            if r != 0:
                return 0 / 0
            cmd = kmc_dump_command + ' -cs' + str(cs_value) + ' '\
                  + os.path.join(work_dir, flist[i][0] + '_sort') + ' '\
                  + os.path.join(output_dir, flist[i][0] + '.txt') + ' '\
                  + os.path.join(output_dir, flist[i][0] + '_beacon.txt')
            r = subprocess.call(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
            if r != 0:
                return 0 / 0
            os.remove(os.path.join(work_dir, flist[i][0] + '.kmc_pre'))
            os.remove(os.path.join(work_dir, flist[i][0] + '.kmc_suf'))
            os.remove(os.path.join(work_dir, flist[i][0] + '_sort.kmc_pre'))
            os.remove(os.path.join(work_dir, flist[i][0] + '_sort.kmc_suf'))
            self.loginfo = 'Number ' + str(countings) + '/' +  str(len(flist))

        return 0
