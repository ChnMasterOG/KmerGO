# coding = utf-8
# author: QiChen
# version: v1.5.0
# modification date: 2020/5/22

import sys, os, shutil, csv
if hasattr(sys, 'frozen'):
    os.environ['PATH'] = sys._MEIPASS + ";" + os.environ['PATH']
from collections import Counter
import platform, ctypes
import multiprocessing
from PyQt5 import QtWidgets, QtCore, QtGui
from qt import MainWindow
from lib import kmc_read, kmer_matrix, kmer_features, sequence_assembly
from lib import projectlist_file as plf

tool_name = 'KmerGO'
project_version = 'V1.5.0'
system_platform = platform.system()

class myWindow(QtWidgets.QMainWindow, MainWindow.Ui_MainWindow):
    def __init__(self):
        # super init
        super(myWindow,self).__init__()
        # class variable
        self.projectfile = plf.ProjectList()
        self.project_dir = ''
        self.GO_flag = False
        self.mode = 1                   # 0-OneClick 1-StepByStep
        self.static_mode = self.mode    # 0-OneClick 1-StepByStep
        # window
        self.new_window = MainWindow.Ui_MainWindow()
        self.new_window.setupUi(self)
        self.new_window.OneClickRunningButton.clicked.connect(self.OneClickRunningButton_Clicked)
        self.new_window.StepByStepRunningButton.clicked.connect(self.StepByStepRunningButton_Clicked)
        self.new_window.Open_Samples_FASTAQ_Button.clicked.connect(self.Open_Samples_FASTAQ_Button_Clicked)
        self.new_window.KMC_Result_Path_Button.clicked.connect(self.KMC_Result_Path_Button_Clicked)
        self.new_window.GM_Result_Path_Button.clicked.connect(self.GM_Result_Path_Button_Clicked)
        self.new_window.GF_Result_Path_Button.clicked.connect(self.GF_Result_Path_Button_Clicked)
        self.new_window.Trait_Info_Path_Button.clicked.connect(self.Trait_Info_Path_Button_Clicked)
        self.new_window.KMC_GO_Button.clicked.connect(self.KMC_GO_Button_Clicked)
        self.new_window.GM_GO_Button.clicked.connect(self.GM_GO_Button_Clicked)
        self.new_window.GF_GO_Button.clicked.connect(self.GF_GO_Button_Clicked)
        self.new_window.KA_GO_Button.clicked.connect(self.KA_GO_Button_Clicked)
        self.new_window.K_Value_Edit.textChanged.connect(self.K_Value_Edit_TextChange)
        self.new_window.CI_Value_Edit.textChanged.connect(self.CI_Value_Edit_TextChange)
        self.new_window.CS_Value_Edit.textChanged.connect(self.CS_Value_Edit_TextChange)
        self.new_window.Process_Number_Edit.textChanged.connect(self.Process_Number_Edit_TextChange)
        self.new_window.ASS_l_Value_Edit.textChanged.connect(self.ASS_l_Value_Edit_TextChange)
        self.new_window.P_Value_Edit.textChanged.connect(self.P_Value_Edit_TextChange)
        self.new_window.ASS_n_Value_Edit.textChanged.connect(self.ASS_n_Value_Edit_TextChange)
        self.new_window.Catagory_RadioButton.clicked.connect(self.Catagory_RadioButton_Clicked)
        self.new_window.Continuous_RadioButton.clicked.connect(self.Continuous_RadioButton_Clicked)
        self.new_window.K_Value_Edit.setToolTip('K-mer length\n(K from 14 to 256; default: 40)')
        self.new_window.K_Value_Label.setToolTip('K-mer length\n(K from 14 to 256; default: 40)')
        self.new_window.CiValue_Label.setToolTip('minimal K-mer occurring times\n(default: 2)')
        self.new_window.CsValue_Label.setToolTip('maximal K-mer occurring times\n(default: 65535)')
        self.new_window.CI_Value_Edit.setToolTip('minimal K-mer occurring times\n(default: 2)')
        self.new_window.CS_Value_Edit.setToolTip('maximal K-mer occurring times\n(default: 65535)')
        self.new_window.Process_Number_Edit.setToolTip('number of processes\n(default: 24)')
        self.new_window.Process_Number_Label.setToolTip('number of processes\n(default: 24)')
        self.new_window.ASS_l_Value_Edit.setToolTip('logical features ASS value\nASS=(TP/P+TN/N)/2\n(default: 0.8)')
        self.new_window.ASS_l_Value_Label.setToolTip('logical features ASS value\nASS=(TP/P+TN/N)/2\n(default: 0.8)')
        self.new_window.P_Value_Edit.setToolTip('numeric features rank sum test p threshold value\n(default: 0.01)')
        self.new_window.P_Value_Label.setToolTip('numeric features rank sum test p threshold value\n(default: 0.01)')
        self.new_window.ASS_n_Value_Edit.setToolTip('numeric features logistic regression ASS value\n(default: 0.8)')
        self.new_window.ASS_n_Value_Label.setToolTip('numeric features logistic regression ASS value\n(default: 0.8)')
        self.setWindowTitle(tool_name + ' ' + project_version)
        self.setFixedSize(self.width(), self.height())
        # initialization
        self.StepByStepRunningButton_Clicked()

    def closeEvent(self, event):
        sys.exit(app.exec_())

    def get_free_space_b(self, folder):
        """ Return folder/drive free space (in bytes) or -1 when this platform is not supported
        """
        if platform.system() == 'Windows':
            free_bytes = ctypes.c_ulonglong(0)
            ctypes.windll.kernel32.GetDiskFreeSpaceExW(ctypes.c_wchar_p(folder), None, None, ctypes.pointer(free_bytes))
            return free_bytes.value
        elif platform.system() == 'Linux':
            st = os.statvfs(folder)
            return st.f_bavail * st.f_frsize
        else:
            return -1

    def get_space_occupation(self, path):
        filesize = 0
        for dirpath, dirname, filename in os.walk(path):
            for ii in filename:
                filesize += os.path.getsize(os.path.join(dirpath,ii))
        return filesize

    def ReadConfiguration(self):
        if self.projectfile.KMC_OK == False:
            self.new_window.KMC_Project_Status_Label.setText('Status: Ready')
            self.new_window.KMC_Project_Status_Label.setStyleSheet('color:red')
        else:
            self.new_window.KMC_Project_Status_Label.setText('Status: Complete')
            self.new_window.KMC_Project_Status_Label.setStyleSheet('color:green')
        if self.projectfile.GM_OK == False:
            self.new_window.GM_Project_Status_Label.setText('Status: Ready')
            self.new_window.GM_Project_Status_Label.setStyleSheet('color:red')
        else:
            self.new_window.GM_Project_Status_Label.setText('Status: Complete')
            self.new_window.GM_Project_Status_Label.setStyleSheet('color:green')
        if self.projectfile.GF_OK == False:
            self.new_window.GF_Project_Status_Label.setText('Status: Ready')
            self.new_window.GF_Project_Status_Label.setStyleSheet('color:red')
        else:
            self.new_window.GF_Project_Status_Label.setText('Status: Complete')
            self.new_window.GF_Project_Status_Label.setStyleSheet('color:green')
        if self.projectfile.KA_OK == False:
            self.new_window.KA_Project_Status_Label.setText('Status: Ready')
            self.new_window.KA_Project_Status_Label.setStyleSheet('color:red')
        else:
            self.new_window.KA_Project_Status_Label.setText('Status: Complete')
            self.new_window.KA_Project_Status_Label.setStyleSheet('color:green')
        self.new_window.Samples_FASTAQ_Path_Edit.setText(self.projectfile.FASTAQ_path)
        self.new_window.K_Value_Edit.setText(str(self.projectfile.K_value))
        self.new_window.CI_Value_Edit.setText(str(self.projectfile.Ci_value))
        self.new_window.CS_Value_Edit.setText(str(self.projectfile.Cs_value))
        self.new_window.KMC_Result_Path_Edit.setText(self.projectfile.KMC_path)
        self.new_window.Process_Number_Edit.setText(str(self.projectfile.Process_value))
        self.new_window.GM_Result_Path_Edit.setText(self.projectfile.GM_path)
        self.new_window.P_Value_Edit.setText(str(self.projectfile.P_value))
        if self.new_window.Catagory_RadioButton.isChecked() is True:
            self.new_window.ASS_n_Value_Edit.setText(str(self.projectfile.ASS_n_value))
        else:
            self.new_window.ASS_n_Value_Edit.setText(str(self.projectfile.Corr_value))
        self.new_window.GF_Result_Path_Edit.setText(self.projectfile.GF_path)
        self.new_window.ASS_l_Value_Edit.setText(str(self.projectfile.ASS_l_value))

    def OneClickRunningButton_Clicked(self):
        getdir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select a folder for one-click runing.', './')
        if getdir == '': return
        self.static_mode = self.mode = 0
        self.new_window.OneClickRunningButton.setText('*One Click Running')
        self.new_window.StepByStepRunningButton.setText('Step By Step Running')
        self.new_window.KA_GO_Button.setText('One-Click Start')
        self.new_window.KA_GO_Button.setGeometry(QtCore.QRect(170, 75, 121, 31))
        self.new_window.KA_GO_Button.setFont(QtGui.QFont("Arial", 11, QtGui.QFont.Bold))
        self.projectfile.KMC_OK = False
        self.projectfile.GM_OK = False
        self.projectfile.GF_OK = False
        self.projectfile.KA_OK = False
        self.new_window.KMC_Result_Path_Button.setEnabled(False)
        self.new_window.GM_Result_Path_Button.setEnabled(False)
        self.new_window.GF_Result_Path_Button.setEnabled(False)
        self.new_window.KMC_GO_Button.setEnabled(False)
        self.new_window.GM_GO_Button.setEnabled(False)
        self.new_window.GF_GO_Button.setEnabled(False)
        self.project_dir = getdir
        if system_platform == 'Windows':
            self.project_dir = self.project_dir.replace('/', '\\')
        elif system_platform == 'Linux':
            self.project_dir = self.project_dir.replace('\\', '/')
        self.setWindowTitle(tool_name + ' ' + project_version + ' Workpath: ' + self.project_dir)
        try:
            os.mkdir(os.path.join(self.project_dir, 'kmer_countings'))
        except:
            pass
        try:
            os.mkdir(os.path.join(self.project_dir, 'kmer_matrix'))
        except:
            pass
        try:
            os.mkdir(os.path.join(self.project_dir, 'kmer_features'))
        except:
            pass
        try:
            os.mkdir(os.path.join(self.project_dir, 'contig_result'))
        except:
            pass
        self.projectfile.KMC_path = os.path.join(self.project_dir, 'kmer_countings')
        self.projectfile.GM_path = os.path.join(self.project_dir, 'kmer_matrix')
        self.projectfile.GF_path = os.path.join(self.project_dir, 'kmer_features')
        self.ReadConfiguration()

    def StepByStepRunningButton_Clicked(self):
        self.static_mode = self.mode = 1
        self.new_window.OneClickRunningButton.setText('One Click Running')
        self.new_window.StepByStepRunningButton.setText('*Step By Step Running')
        self.new_window.KA_GO_Button.setText('Start')
        self.new_window.KA_GO_Button.setGeometry(QtCore.QRect(350, 80, 101, 21))
        self.new_window.KA_GO_Button.setFont(QtGui.QFont("Arial", 10, QtGui.QFont.Normal))
        self.projectfile.KMC_OK = False
        self.projectfile.GM_OK = False
        self.projectfile.GF_OK = False
        self.projectfile.KA_OK = False
        self.new_window.KMC_Result_Path_Button.setEnabled(True)
        self.new_window.GM_Result_Path_Button.setEnabled(True)
        self.new_window.GF_Result_Path_Button.setEnabled(True)
        self.new_window.KMC_GO_Button.setEnabled(True)
        self.new_window.GM_GO_Button.setEnabled(True)
        self.new_window.GF_GO_Button.setEnabled(True)
        self.setWindowTitle(tool_name + ' ' + project_version + ' Workpath: ./')
        self.ReadConfiguration()

    def Open_Samples_FASTAQ_Button_Clicked(self):
        getdir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select the folder of samples sequencing files',
                                                            self.projectfile.FASTAQ_path)
        if getdir == '': return
        if system_platform == 'Windows':
            getdir = getdir.replace('/', '\\')
        elif system_platform == 'Linux':
            getdir = getdir.replace('\\', '/')
        flist = os.listdir(getdir)
        for seq in list(flist):
            suffix = seq[seq.find('.'):]
            if suffix.lower() not in kmc_read.FASTA_suffix and suffix.lower() not in kmc_read.FASTQ_suffix:
                flist.remove(seq)
        self.projectfile.FASTAQ_path = getdir
        self.new_window.Samples_FASTAQ_Path_Edit.setText(self.projectfile.FASTAQ_path)

    def KMC_Result_Path_Button_Clicked(self):
        getdir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select the folder of k-mer counting files',
                                                            self.projectfile.KMC_path)
        if getdir == '': return
        if system_platform == 'Windows':
            getdir = getdir.replace('/', '\\')
        elif system_platform == 'Linux':
            getdir = getdir.replace('\\', '/')
        self.projectfile.KMC_path = getdir
        self.new_window.KMC_Result_Path_Edit.setText(self.projectfile.KMC_path)

    def GM_Result_Path_Button_Clicked(self):
        getdir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select the folder of matrix files',
                                                            self.projectfile.GM_path)
        if getdir == '': return
        if system_platform == 'Windows':
            getdir = getdir.replace('/', '\\')
        elif system_platform == 'Linux':
            getdir = getdir.replace('\\', '/')
        self.projectfile.GM_path = getdir
        self.new_window.GM_Result_Path_Edit.setText(self.projectfile.GM_path)

    def GF_Result_Path_Button_Clicked(self):
        getdir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select the folder of feature files',
                                                            self.projectfile.GF_path)
        if getdir == '': return
        if system_platform == 'Windows':
            getdir = getdir.replace('/', '\\')
        elif system_platform == 'Linux':
            getdir = getdir.replace('\\', '/')
        self.projectfile.GF_path = getdir
        self.new_window.GF_Result_Path_Edit.setText(self.projectfile.GF_path)

    def Trait_Info_Path_Button_Clicked(self):
        getdir = QtWidgets.QFileDialog.getOpenFileName(self, 'Select the csv file of trait information',
                                                       self.projectfile.GF_path, 'csv files(*.csv)')
        getdir = getdir[0]
        if getdir == '': return
        if system_platform == 'Windows':
            getdir = getdir.replace('/', '\\')
        elif system_platform == 'Linux':
            getdir = getdir.replace('\\', '/')
        try:
            csv_f = open(getdir, 'r')
            reader = csv.reader(csv_f)
            for row in reader:
                if row[1] == 'trait':
                    continue
                else:
                    self.projectfile.TI_dic[row[0]] = row[1]
        except:
            QtWidgets.QMessageBox.critical(self, 'Error', 'Can not read the CSV file!', QtWidgets.QMessageBox.Yes)
            return
        if len(self.projectfile.TI_dic) == 0:
            QtWidgets.QMessageBox.critical(self, 'Error', 'Can not read the CSV file!', QtWidgets.QMessageBox.Yes)
            return
        self.projectfile.TI_path = getdir
        self.new_window.Trait_Info_Path_Edit.setText(self.projectfile.TI_path)

    def K_Value_Edit_TextChange(self):
        try: self.projectfile.K_value = int(self.new_window.K_Value_Edit.text())
        except: pass
        self.new_window.K_Value_Edit.setText(str(self.projectfile.K_value))

    def CI_Value_Edit_TextChange(self):
        try: self.projectfile.Ci_value = int(self.new_window.CI_Value_Edit.text())
        except: pass
        self.new_window.CI_Value_Edit.setText(str(self.projectfile.Ci_value))

    def CS_Value_Edit_TextChange(self):
        try: self.projectfile.Cs_value = int(self.new_window.CS_Value_Edit.text())
        except: pass
        self.new_window.CS_Value_Edit.setText(str(self.projectfile.Cs_value))

    def Process_Number_Edit_TextChange(self):
        try: self.projectfile.Process_value = int(self.new_window.Process_Number_Edit.text())
        except: pass
        self.new_window.Process_Number_Edit.setText(str(self.projectfile.Process_value))

    def ASS_l_Value_Edit_TextChange(self):
        try: self.projectfile.ASS_l_value = float(self.new_window.ASS_l_Value_Edit.text())
        except: pass
        self.new_window.ASS_l_Value_Edit.setText(str(self.projectfile.ASS_l_value))

    def P_Value_Edit_TextChange(self):
        try: self.projectfile.P_value = float(self.new_window.P_Value_Edit.text())
        except: pass
        self.new_window.P_Value_Edit.setText(str(self.projectfile.P_value))

    def ASS_n_Value_Edit_TextChange(self):
        if self.new_window.Catagory_RadioButton.isChecked() is True:
            try: self.projectfile.ASS_n_value = float(self.new_window.ASS_n_Value_Edit.text())
            except: pass
            self.new_window.ASS_n_Value_Edit.setText(str(self.projectfile.ASS_n_value))
        else:
            try: self.projectfile.Corr_value = float(self.new_window.ASS_n_Value_Edit.text())
            except: pass
            self.new_window.ASS_n_Value_Edit.setText(str(self.projectfile.Corr_value))

    def KMC_GO_Button_Clicked(self):
        if self.GO_flag == True:
            QtWidgets.QMessageBox.critical(self, 'Error', 'A task is running!', QtWidgets.QMessageBox.Yes)
            return
        self.GO_flag = True
        self.new_window.KMC_GO_Button.setEnabled(False)
        self.new_window.Open_Samples_FASTAQ_Button.setEnabled(False)
        self.new_window.KMC_Result_Path_Button.setEnabled(False)
        self.new_window.K_Value_Edit.setEnabled(False)
        self.new_window.CI_Value_Edit.setEnabled(False)
        self.new_window.CS_Value_Edit.setEnabled(False)
        self.new_window.OneClickRunningButton.setEnabled(False)
        self.new_window.StepByStepRunningButton.setEnabled(False)
        self.kmc_thread = kmc_read.KMC_Thread((self.projectfile.K_value, self.projectfile.Ci_value,
                                               self.projectfile.Cs_value, self.projectfile.FASTAQ_path,
                                               self.projectfile.KMC_path))
        self.kmc_thread.start()
        self.kmc_timer = QtCore.QTimer(self)
        self.kmc_timer.timeout.connect(self.KMC_Timer_Show)
        self.kmc_timer.start(100)

    def KMC_Timer_Show(self):
        if self.kmc_thread.status <= 0:
            self.GO_flag = False
            self.new_window.Open_Samples_FASTAQ_Button.setEnabled(True)
            if self.mode == 1:
                self.new_window.KMC_Result_Path_Button.setEnabled(True)
                self.new_window.KMC_GO_Button.setEnabled(True)
            elif self.kmc_thread.status < 0:
                self.new_window.KA_GO_Button.setEnabled(True)
                self.new_window.Catagory_RadioButton.setEnabled(True)
                self.new_window.Continuous_RadioButton.setEnabled(True)
                self.new_window.Trait_Info_Path_Button.setEnabled(True)
            self.new_window.K_Value_Edit.setEnabled(True)
            self.new_window.CI_Value_Edit.setEnabled(True)
            self.new_window.CS_Value_Edit.setEnabled(True)
            self.new_window.OneClickRunningButton.setEnabled(True)
            self.new_window.StepByStepRunningButton.setEnabled(True)
            if self.kmc_thread.status == 0:
                self.projectfile.KMC_OK = True
                self.new_window.KMC_Project_Status_Label.setText('Status: Complete')
                self.new_window.KMC_Project_Status_Label.setStyleSheet('color:green')
            else:
                self.new_window.KMC_Project_Status_Label.setText('Status: ' + self.kmc_thread.loginfo)
                self.new_window.KMC_Project_Status_Label.setStyleSheet('color:red')
            self.kmc_timer.stop()
            # if one-click
            if self.kmc_thread.status == 0 and self.mode == 0:
                del self.kmc_thread  # free up memery
                self.GM_GO_Button_Clicked()
            else:
                del self.kmc_thread  # free up memery
        else:
            self.new_window.KMC_Project_Status_Label.setText('Status: ' + self.kmc_thread.loginfo)
            self.new_window.KMC_Project_Status_Label.setStyleSheet('color:blue')

    def GM_GO_Button_Clicked(self):
        if self.GO_flag == True:
            QtWidgets.QMessageBox.critical(self, 'Error', 'A tast is running!', QtWidgets.QMessageBox.Yes)
            return
        # check free space
        kmc_file_size = self.get_space_occupation(self.projectfile.KMC_path)
        free_space_size = self.get_free_space_b(self.projectfile.KMC_path)
        Group_counting = self.projectfile.GroupA_Number + self.projectfile.GroupB_Number
        file_size_predicted = kmc_file_size / (self.projectfile.K_value + len(str(self.projectfile.Cs_value)) + 2) / (0.2 * Group_counting) \
                              * (self.projectfile.K_value + 1.8 * Group_counting + 0.2 * Group_counting * len(str(self.projectfile.Cs_value)) - 1)
        if file_size_predicted > free_space_size:
            buttonReply = QtWidgets.QMessageBox.warning(self, 'Space is not enough', 
                                                        'You only have %.4fM disk space left, but maybe result files will over %.4fM.\nDo you want to continue?' 
                                                        % (free_space_size / 1024 / 1024, file_size_predicted / 1024 / 1024), 
                                                        QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
            if buttonReply != True:
                return
        # initialization folders
        try:
            os.remove('temp')
        except:
            pass
        self.GO_flag = True
        self.new_window.KMC_Result_Path_Button.setEnabled(False)
        self.new_window.GM_Result_Path_Button.setEnabled(False)
        self.new_window.Process_Number_Edit.setEnabled(False)
        self.new_window.GM_GO_Button.setEnabled(False)
        self.new_window.OneClickRunningButton.setEnabled(False)
        self.new_window.StepByStepRunningButton.setEnabled(False)
        if self.static_mode == 0 and self.new_window.Catagory_RadioButton.isChecked() is True:
            self.gm_thread = kmer_matrix.GM_Thread((self.projectfile.KMC_path, self.projectfile.GM_path,
                                                    self.projectfile.Process_value, self.projectfile.GroupA_Name,
                                                    self.projectfile.GroupA_Number, self.projectfile.TI_dic))
        elif self.static_mode == 0 and self.new_window.Catagory_RadioButton.isChecked() is False:
            self.gm_thread = kmer_matrix.GM_Thread((self.projectfile.KMC_path, self.projectfile.GM_path,
                                                    self.projectfile.Process_value, '',
                                                    len(self.projectfile.TI_dic), self.projectfile.TI_dic))
        else:
            self.gm_thread = kmer_matrix.GM_Thread((self.projectfile.KMC_path, self.projectfile.GM_path,
                                                   self.projectfile.Process_value, None, None, None))
        self.gm_thread.start()
        self.gm_timer = QtCore.QTimer(self)
        self.gm_timer.timeout.connect(self.GM_Timer_Show)
        self.gm_timer.start(1000)

    def GM_Timer_Show(self):
        if self.gm_thread.status <= 0:
            self.GO_flag = False
            if self.mode == 1:
                self.new_window.KMC_Result_Path_Button.setEnabled(True)
                self.new_window.GM_Result_Path_Button.setEnabled(True)
                self.new_window.GM_GO_Button.setEnabled(True)
            elif self.gm_thread.status < 0:
                self.new_window.KA_GO_Button.setEnabled(True)
                self.new_window.Catagory_RadioButton.setEnabled(True)
                self.new_window.Continuous_RadioButton.setEnabled(True)
                self.new_window.Trait_Info_Path_Button.setEnabled(True)
            self.new_window.Process_Number_Edit.setEnabled(True)
            self.new_window.OneClickRunningButton.setEnabled(True)
            self.new_window.StepByStepRunningButton.setEnabled(True)
            try:
                shutil.rmtree('temp')
                for job in self.gm_thread.jobs:
                    job.terminate()
            except:
                pass
            if self.gm_thread.status == 0:
                self.projectfile.GM_OK = True
                self.new_window.GM_Project_Status_Label.setText('Status: Complete')
                self.new_window.GM_Project_Status_Label.setStyleSheet('color:green')
                self.new_window.GM_Project_Status_Label.setToolTip('')
            else:
                self.new_window.GM_Project_Status_Label.setText('Status: ' + self.gm_thread.loginfo)
                self.new_window.GM_Project_Status_Label.setStyleSheet('color:red')
            self.gm_timer.stop()
            # if one-click
            if self.gm_thread.status == 0 and self.mode == 0:
                del self.gm_thread  # free up memery
                self.GF_GO_Button_Clicked()
            else:
                del self.gm_thread  # free up memery
        else:
            if self.gm_thread.status == 2:
                error_status = self.gm_thread.detective_error()
                if error_status == 0:
                    self.new_window.GM_Project_Status_Label.setText('Status: Move the mouse here to see the details')
                    temp_progress = self.gm_thread.detective_progress()
                    tipstr = 'Total: %.2f%%\n' % (100 * sum(temp_progress) /
                                                  sum(self.gm_thread.filesize))
                    for k in range(self.gm_thread.process_number):
                        if self.gm_thread.block_size[k] == 0:
                            tipstr += 'Process ' + str(k + 1) + ': 100%\n'
                        else:
                            tipstr += 'Process ' + str(k + 1) + ': %.2f%%\n' % \
                                      (100 * temp_progress[k] / self.gm_thread.block_size[k])
                    self.new_window.GM_Project_Status_Label.setToolTip(tipstr[:-1])
                    self.new_window.GM_Project_Status_Label.setStyleSheet('color:blue')
                else:
                    if error_status == -1:
                        self.gm_thread.loginfo = 'Missing some KMC result files.'
                    elif error_status == -2:
                        self.gm_thread.loginfo = 'Missing the result folder.'
                    elif error_status == -3:
                        self.gm_thread.loginfo = 'Something is wrong when files\' pointer move.'
                    elif error_status == -10:
                        self.gm_thread.loginfo = 'Parameters error? Processes break down.'
                    self.gm_thread.status = -9

    def GF_GO_Button_Clicked(self):
        if self.GO_flag == True:
            QtWidgets.QMessageBox.critical(self, 'Error', 'A task is running!', QtWidgets.QMessageBox.Yes)
            return

        if self.mode != 0 and self.Check_csv_validity() == False:
            QtWidgets.QMessageBox.critical(self, 'Error', 'CSV format error', QtWidgets.QMessageBox.Yes)
            return

        # check free space
        gm_file_size = self.get_space_occupation(self.projectfile.GM_path)
        free_space_size = self.get_free_space_b(self.projectfile.GM_path)
        file_size_predicted = gm_file_size
        if file_size_predicted > free_space_size:
            buttonReply = QtWidgets.QMessageBox.warning(self, 'Space is not enough', 
                                                        'You only have %.4fM disk space left, but maybe result files will over %.4fM.\nDo you want to continue?'
                                                        % (free_space_size / 1024 / 1024, file_size_predicted / 1024 / 1024), 
                                                        QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
            if buttonReply != True:
                return
        # initialization folders
        try:
            os.remove('temp')
        except:
            pass
        self.GO_flag = True
        self.new_window.GF_GO_Button.setEnabled(False)
        self.new_window.GM_Result_Path_Button.setEnabled(False)
        self.new_window.GF_Result_Path_Button.setEnabled(False)
        self.new_window.ASS_l_Value_Edit.setEnabled(False)
        self.new_window.P_Value_Edit.setEnabled(False)
        self.new_window.ASS_n_Value_Edit.setEnabled(False)
        self.new_window.OneClickRunningButton.setEnabled(False)
        self.new_window.StepByStepRunningButton.setEnabled(False)
        self.new_window.Catagory_RadioButton.setEnabled(False)
        self.new_window.Continuous_RadioButton.setEnabled(False)
        self.new_window.Trait_Info_Path_Button.setEnabled(False)
        self.gf_thread = kmer_features.GF_Thread((self.projectfile.GM_path, self.projectfile.GF_path,
                                                  self.projectfile.ASS_l_value, self.projectfile.P_value,
                                                  self.projectfile.ASS_n_value, self.projectfile.GroupA_Number,
                                                  self.projectfile.GroupB_Number, self.projectfile.GroupA_Name,
                                                  self.projectfile.GroupB_Name, self.projectfile.TI_dic,
                                                  self.projectfile.Corr_value,
                                                  self.new_window.Catagory_RadioButton.isChecked()))
        self.gf_thread.start()
        self.gf_timer = QtCore.QTimer(self)
        self.gf_timer.timeout.connect(self.GF_Timer_Show)
        self.gf_timer.start(1000)

    def GF_Timer_Show(self):
        if self.gf_thread.status <= 0:
            self.GO_flag = False
            if self.mode == 1:
                self.new_window.GM_Result_Path_Button.setEnabled(True)
                self.new_window.GF_Result_Path_Button.setEnabled(True)
                self.new_window.GF_GO_Button.setEnabled(True)
            elif self.gf_thread.status < 0:
                self.new_window.KA_GO_Button.setEnabled(True)
                self.new_window.Catagory_RadioButton.setEnabled(True)
                self.new_window.Continuous_RadioButton.setEnabled(True)
                self.new_window.Trait_Info_Path_Button.setEnabled(True)
            self.new_window.ASS_l_Value_Edit.setEnabled(True)
            self.new_window.P_Value_Edit.setEnabled(True)
            self.new_window.ASS_n_Value_Edit.setEnabled(True)
            self.new_window.OneClickRunningButton.setEnabled(True)
            self.new_window.StepByStepRunningButton.setEnabled(True)
            if self.mode == 1:
                self.new_window.Catagory_RadioButton.setEnabled(True)
                self.new_window.Continuous_RadioButton.setEnabled(True)
                self.new_window.Trait_Info_Path_Button.setEnabled(True)
            try:
                shutil.rmtree('temp')
                for job in self.gf_thread.jobs:
                    job.terminate()
            except:
                pass
            if self.gf_thread.status == 0:
                self.projectfile.GF_OK = True
                self.new_window.GF_Project_Status_Label.setText('Status: Complete')
                self.new_window.GF_Project_Status_Label.setStyleSheet('color:green')
                self.new_window.GF_Project_Status_Label.setToolTip('')
            else:
                self.new_window.GF_Project_Status_Label.setText('Status: ' + self.gf_thread.loginfo)
                self.new_window.GF_Project_Status_Label.setStyleSheet('color:red')
            self.gf_timer.stop()
            # if one-click
            if self.gf_thread.status == 0 and self.mode == 0:
                self.mode = 1
                del self.gf_thread  # free up memery
                self.KA_GO_Button_Clicked()
            else:
                del self.gf_thread  # free up memery
        else:
            if self.gf_thread.status == 2:
                error_status = self.gf_thread.detective_error()
                if error_status == 0:
                    self.new_window.GF_Project_Status_Label.setText('Status: Move the mouse here to see the details')
                    temp_progress = self.gf_thread.detective_progress()
                    tipstr = 'Total: %.2f%%\n' % (100 * sum(temp_progress) /
                                                  sum(self.gf_thread.filesize))
                    for k in range(self.gf_thread.files_number):
                        tipstr += 'Process ' + str(k + 1) + ': %.2f%%\n' % \
                                  (100 * temp_progress[k] / self.gf_thread.filesize[k])
                    self.new_window.GF_Project_Status_Label.setToolTip(tipstr[:-1])
                    self.new_window.GF_Project_Status_Label.setStyleSheet('color:blue')
                else:
                    if error_status == -1:
                        self.gf_thread.loginfo = 'Missing matrix files.'
                    elif error_status == -2:
                        self.gf_thread.loginfo = 'Missing the result folder.'
                    self.gf_thread.status = -9

    def KA_GO_Button_Clicked(self):
        if self.GO_flag == True:
            QtWidgets.QMessageBox.critical(self, 'Error', 'A task is running!', QtWidgets.QMessageBox.Yes)
            return
        # one-click running or step-by-step running
        if self.mode == 0:
            if self.projectfile.TI_path == '':
                QtWidgets.QMessageBox.critical(self, 'Error', 'No CSV files selected!', QtWidgets.QMessageBox.Yes)
                return
            if self.Check_csv_validity() == False:
                QtWidgets.QMessageBox.critical(self, 'Error', 'CSV format error', QtWidgets.QMessageBox.Yes)
                return
            self.new_window.KA_GO_Button.setEnabled(False)
            self.new_window.Catagory_RadioButton.setEnabled(False)
            self.new_window.Continuous_RadioButton.setEnabled(False)
            self.new_window.Trait_Info_Path_Button.setEnabled(False)
            self.KMC_GO_Button_Clicked()
            return
        # initialization folders
        try:
            os.remove('temp')
        except:
            pass
        self.GO_flag = True
        self.new_window.KA_GO_Button.setEnabled(False)
        self.new_window.GF_Result_Path_Button.setEnabled(False)
        self.new_window.OneClickRunningButton.setEnabled(False)
        self.new_window.StepByStepRunningButton.setEnabled(False)
        if not(os.path.exists(os.path.join(self.project_dir, 'contig_result'))):
            os.mkdir(os.path.join(self.project_dir, 'contig_result'))
        self.ka_thread = sequence_assembly.KA_Thread((self.projectfile.GF_path,
                                                      os.path.join(self.project_dir, 'contig_result')))
        self.ka_thread.start()
        self.ka_timer = QtCore.QTimer(self)
        self.ka_timer.timeout.connect(self.KA_Timer_Show)
        self.ka_timer.start(100)

    def KA_Timer_Show(self):
        if self.ka_thread.status <= 0:
            self.GO_flag = False
            self.mode = self.static_mode
            self.new_window.KA_GO_Button.setEnabled(True)
            if self.mode == 1:
                self.new_window.GF_Result_Path_Button.setEnabled(True)
            self.new_window.Catagory_RadioButton.setEnabled(True)
            self.new_window.Continuous_RadioButton.setEnabled(True)
            self.new_window.Trait_Info_Path_Button.setEnabled(True)
            self.new_window.OneClickRunningButton.setEnabled(True)
            self.new_window.StepByStepRunningButton.setEnabled(True)
            if self.ka_thread.status == 0:
                self.projectfile.KA_OK = True
                self.new_window.KA_Project_Status_Label.setText('Status: Complete')
                self.new_window.KA_Project_Status_Label.setStyleSheet('color:green')
                if self.mode == 1:
                    QtWidgets.QMessageBox.information(self, 'Success', 'Result files are stored in: ' +
                                                      os.path.join(self.project_dir, 'contig_result'),
                                                      QtWidgets.QMessageBox.Yes)
                else:
                    QtWidgets.QMessageBox.information(self, 'Success', 'Done!\nResult files are stored in: ' +
                                                      os.path.join(self.project_dir, 'contig_result'),
                                                      QtWidgets.QMessageBox.Yes)
            else:
                self.new_window.KA_Project_Status_Label.setText('Status: ' + self.ka_thread.loginfo)
                self.new_window.KA_Project_Status_Label.setStyleSheet('color:red')
            del self.ka_thread     # free up memery
            self.ka_timer.stop()
        else:
            self.new_window.KA_Project_Status_Label.setText('Status: ' + self.ka_thread.loginfo)
            self.new_window.KA_Project_Status_Label.setStyleSheet('color:blue')

    def Catagory_RadioButton_Clicked(self):
        self.new_window.ASS_l_Value_Label.setVisible(True)
        self.new_window.ASS_l_Value_Edit.setVisible(True)
        self.new_window.ASS_n_Value_Label.setText('<html><head/><body><p align="center">ASS-n =</p></body></html>')
        self.new_window.ASS_n_Value_Edit.setText(str(self.projectfile.ASS_n_value))
        self.new_window.P_Value_Edit.setToolTip('numeric features rank sum test p threshold value\n(default: 0.01)')
        self.new_window.P_Value_Label.setToolTip('numeric features rank sum test p threshold value\n(default: 0.01)')
        self.new_window.ASS_n_Value_Edit.setToolTip('numeric features logistic regression ASS value\n(default: 0.8)')
        self.new_window.ASS_n_Value_Label.setToolTip('numeric features logistic regression ASS value\n(default: 0.8)')

    def Continuous_RadioButton_Clicked(self):
        self.new_window.ASS_l_Value_Label.setVisible(False)
        self.new_window.ASS_l_Value_Edit.setVisible(False)
        self.new_window.ASS_n_Value_Label.setText('<html><head/><body><p align="center">ρ =</p></body></html>')
        self.new_window.ASS_n_Value_Edit.setText(str(self.projectfile.Corr_value))
        self.new_window.P_Value_Edit.setToolTip('logical features rank sum test p threshold value\n(default: 0.01)')
        self.new_window.P_Value_Label.setToolTip('logical features rank sum test p threshold value\n(default: 0.01)')
        self.new_window.ASS_n_Value_Edit.setToolTip('numeric features coefficient of association ρ threshold value\n(default: 0.8)')
        self.new_window.ASS_n_Value_Label.setToolTip('numeric features coefficient of association ρ threshold value\n(default: 0.8)')

    def Check_csv_validity(self):
        if self.new_window.Catagory_RadioButton.isChecked() is True:
            counter_dic = Counter(self.projectfile.TI_dic.values())
            if len(counter_dic) != 2:
                return False
            else:
                self.projectfile.GroupA_Name = list(counter_dic.keys())[0]
                self.projectfile.GroupA_Number = list(counter_dic.values())[0]
                self.projectfile.GroupB_Name = list(counter_dic.keys())[1]
                self.projectfile.GroupB_Number = list(counter_dic.values())[1]
                self.projectfile.Group_Number = self.projectfile.GroupA_Number + self.projectfile.GroupB_Number
                return True
        else:
            for v in self.projectfile.TI_dic.values():
                try:
                    float(v)
                except:
                    return False
            return True

if __name__ == "__main__":
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    multiprocessing.freeze_support()
    app = QtWidgets.QApplication(sys.argv)
    mywindow = myWindow()
    mywindow.show()
    sys.exit(app.exec_())
