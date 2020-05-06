# coding = utf-8
# author: QiChen
# version: v2.0
# modification date: 2020/4/26

import datetime, os

class ProjectList:
    def __init__(self):
        self.FASTAQ_path = ''
        self.KMC_OK = False
        self.KMC_path = ''
        self.GM_OK = False
        self.GM_path = ''
        self.GF_OK = False
        self.GF_path = ''
        self.TI_path = ''
        self.GF_mode = 0
        self.TI_dic = {}
        self.KA_OK = False
        self.K_value = 40
        self.Ci_value = 2
        self.Cs_value = 65536
        self.Process_value = 24
        self.ASS_l_value = 0.8
        self.P_value = 0.01
        self.ASS_n_value = 0.8
        self.Corr_value = 0.8
        self.GroupA_Number = 1  # Set up automatically according to fasta/q files number
        self.GroupB_Number = 1  # Set up automatically according to fasta/q file number
        self.GroupA_Name = ''
        self.GroupB_Name = ''
        self.Group_Number = 1
        self.Modification_record = ''

    def CreateNewFile(self, project_dir):
        self.__init__()
        self.KMC_path = os.path.join(project_dir, 'kmer_countings')
        self.GM_path = os.path.join(project_dir, 'kmer_matrix')
        self.GF_path = os.path.join(project_dir, 'kmer_features')
        self.WriteFile(project_dir)

    def WriteFile(self, project_dir):
        now = datetime.datetime.now()
        self.Modification_record += now.strftime('%Y-%m-%d %H:%M:%S') + '\n'
        f = open(os.path.join(project_dir, 'ProjectList.list'), 'w')
        f.write('--------------------Project Parameters--------------------\n')
        f.write('FASTAQ_path = ' + self.FASTAQ_path + '\n')
        f.write('KMC_OK = ' + str(self.KMC_OK) + '\n')
        f.write('KMC_path = ' + self.KMC_path + '\n')
        f.write('GM_OK = ' + str(self.GM_OK) + '\n')
        f.write('GM_path = ' + self.GM_path + '\n')
        f.write('GF_OK = ' + str(self.GF_OK) + '\n')
        f.write('GF_path = ' + self.GF_path + '\n')
        f.write('KA_OK = ' + str(self.KA_OK) + '\n')
        f.write('K_value = ' + str(self.K_value) + '\n')
        f.write('Ci_value = ' + str(self.Ci_value) + '\n')
        f.write('Cs_value = ' + str(self.Cs_value) + '\n')
        f.write('Process_value = ' + str(self.Process_value) + '\n')
        f.write('ASS_l_value = ' + str(self.ASS_l_value) + '\n')
        f.write('P_value = ' + str(self.P_value) + '\n')
        f.write('ASS_n_value = ' + str(self.ASS_n_value) + '\n')
        f.write('GroupA_Number = ' + str(self.GroupA_Number) + '\n')
        f.write('GroupB_Number = ' + str(self.GroupB_Number) + '\n')
        f.write('--------------------Modification History--------------------\n')
        f.write(self.Modification_record)
        f.close()

    def ReadFile(self, project_dir):
        f = open(os.path.join(project_dir, 'ProjectList.list'), 'r')
        f.readline()    # skip project parameters
        self.FASTAQ_path = f.readline()[16:-1]
        self.KMC_OK = (lambda x: True if x == 'True' else False)(f.readline()[9:-1])
        self.KMC_path = f.readline()[11:-1]
        self.GM_OK = (lambda x: True if x == 'True' else False)(f.readline()[8:-1])
        self.GM_path = f.readline()[10:-1]
        self.GF_OK = (lambda x: True if x == 'True' else False)(f.readline()[8:-1])
        self.GF_path = f.readline()[10:-1]
        self.KA_OK = (lambda x: True if x == 'True' else False)(f.readline()[8:-1])
        self.K_value = int(f.readline()[10:-1])
        self.Ci_value = int(f.readline()[11:-1])
        self.Cs_value = int(f.readline()[11:-1])
        self.Process_value = int(f.readline()[16:-1])
        self.ASS_l_value = float(f.readline()[14:-1])
        self.P_value = float(f.readline()[10:-1])
        self.ASS_n_value = float(f.readline()[14:-1])
        self.GroupA_Number = int(f.readline()[16:-1])
        self.GroupB_Number = int(f.readline()[16:-1])
        f.readline()    # skip modification history
        self.Modification_record = f.read()
        f.close()
