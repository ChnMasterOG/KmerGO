# coding = utf-8
# author: QiChen
# version: v3.1
# modification date: 2021/9/25

import os, sys, shutil
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
        kmc_input = []
        if len(flist) == 0:
            return -2  # no files in this directory
        for i in list(flist):
            dot_pos = i.find('.')
            suffix = i[dot_pos:]
            is_pairend = i[dot_pos-2:dot_pos] == "_1"
            if suffix.lower() in FASTA_suffix:
                typelist.append(' -fm ')
            elif suffix.lower() in FASTQ_suffix:
                typelist.append(' -fq ')
            else:
                continue
            if is_pairend == False:
                kmc_input.append(i)
            else:
                typelist = typelist[:-1]
                pair_1_name = i[:dot_pos-2] + i[dot_pos:]
                if pair_1_name in kmc_input:
                    kmc_input[kmc_input.index(pair_1_name)] = "@" + i[:dot_pos-2] + ".list"
                    listf = open(os.path.join(work_dir, i[:dot_pos-2] + ".list"), 'w')
                    listf.write(os.path.join(work_dir, pair_1_name) + '\n')
                    listf.write(os.path.join(work_dir, i) + '\n')
                    listf.close()

        countings = 0
        self.loginfo = 'Number 0/' + str(len(kmc_input))
        for i in range(len(kmc_input)):
            if kmc_input[i][0] == '@':
                prechar = '@'
                kmc_input[i] = kmc_input[i][1:]
            else:
                prechar = ''
            kmerdb_name = kmc_input[i][:kmc_input[i].find('.')]
            countings += 1
            cmd = kmc_command + ' -k' + str(k_value) + ' -ci' + str(ci_value) + ' -cs' + str(cs_value)\
                  + typelist[i] + prechar + os.path.join(work_dir, kmc_input[i]) + ' ' + os.path.join(work_dir, kmerdb_name)\
                  + ' ' + work_dir
            r = subprocess.call(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
            if r != 0:
                return 0 / 0
            if k_value >= 14:
                cmd = kmc_tools_command + ' transform ' + os.path.join(work_dir, kmerdb_name)\
                      + ' sort ' + os.path.join(work_dir, kmerdb_name + '_sort')
                r = subprocess.call(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
                if r != 0:
                    return 0 / 0
            else:
                shutil.copyfile(os.path.join(work_dir, kmerdb_name + '.kmc_pre'), os.path.join(work_dir, kmerdb_name + '_sort.kmc_pre'))
                shutil.copyfile(os.path.join(work_dir, kmerdb_name + '.kmc_suf'), os.path.join(work_dir, kmerdb_name + '_sort.kmc_suf'))
            cmd = kmc_dump_command + ' -cs' + str(cs_value) + ' '\
                  + os.path.join(work_dir, kmerdb_name + '_sort') + ' '\
                  + os.path.join(output_dir, kmerdb_name + '.txt') + ' '\
                  + os.path.join(output_dir, kmerdb_name + '_beacon.txt')
            r = subprocess.call(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
            if r != 0:
                return 0 / 0
            os.remove(os.path.join(work_dir, kmerdb_name + '.kmc_pre'))
            os.remove(os.path.join(work_dir, kmerdb_name + '.kmc_suf'))
            os.remove(os.path.join(work_dir, kmerdb_name + '_sort.kmc_pre'))
            os.remove(os.path.join(work_dir, kmerdb_name + '_sort.kmc_suf'))
            self.loginfo = 'Number ' + str(countings) + '/' +  str(len(kmc_input))

            # clear listf
            if prechar == '@':
                os.remove(os.path.join(work_dir, kmc_input[i]))

        return 0
