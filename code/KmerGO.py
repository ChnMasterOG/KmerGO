# coding = utf-8
# author: QiChen
# version: v1.3.0
# modification date: 2019/12/10

import sys, os, shutil
if hasattr(sys, 'frozen'):
    os.environ['PATH'] = sys._MEIPASS + ";" + os.environ['PATH']
import platform, ctypes
from PyQt5 import QtWidgets, QtCore
from qt import MainWindow
from lib import kmc_read, kmer_matrix, kmer_features, sequence_assembly
from lib import projectlist_file as plf

tool_name = 'KmerGO'
project_version = 'V1.3.0'
system_platform = platform.system()

class myWindow(QtWidgets.QMainWindow, MainWindow.Ui_MainWindow):
    def __init__(self):
        # super init
        super(myWindow,self).__init__()
        # class variable
        self.projectfile = plf.ProjectList()
        self.project_dir = ''
        self.GO_flag = False
        # window
        self.new_window = MainWindow.Ui_MainWindow()
        self.new_window.setupUi(self)
        self.new_window.NewProjectButton.clicked.connect(self.NewProject_Clicked)
        self.new_window.OpenProjectButton.clicked.connect(self.OpenProject_Clicked)
        self.new_window.SaveProjectButton.clicked.connect(self.SaveProject_Clicked)
        self.new_window.Open_GroupA_FASTAQ_Button.clicked.connect(self.Open_GroupA_FASTAQ_Button_Clicked)
        self.new_window.Open_GroupB_FASTAQ_Button.clicked.connect(self.Open_GroupB_FASTAQ_Button_Clicked)
        self.new_window.KMC_Result_Path_Button.clicked.connect(self.KMC_Result_Path_Button_Clicked)
        self.new_window.GM_Result_Path_Button.clicked.connect(self.GM_Result_Path_Button_Clicked)
        self.new_window.GF_Result_Path_Button.clicked.connect(self.GF_Result_Path_Button_Clicked)
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
        self.new_window.GroupA_Number_Edit.textChanged.connect(self.GroupA_Number_Edit_TextChange)
        self.new_window.GroupB_Number_Edit.textChanged.connect(self.GroupB_Number_Edit_TextChange)
        self.new_window.K_Value_Edit.setToolTip('K-mer length\n(K from 1 to 256; default: 40)')
        self.new_window.CI_Value_Edit.setToolTip('exclude K-mers occurring less than <value> times\n(default: 2)')
        self.new_window.CS_Value_Edit.setToolTip('maximal value of a counter\n(default: 65535)')
        self.new_window.Process_Number_Edit.setToolTip('number of processes\n(default: 24)')
        self.new_window.ASS_l_Value_Edit.setToolTip('logical features ASS value\nASS=(TP/P+TN/N)/2\n(default: 0.8)')
        self.new_window.P_Value_Edit.setToolTip('numeric features rank sum test p threshold value\n(default: 0.01)')
        self.new_window.ASS_n_Value_Edit.setToolTip('numeric features logistic regression ASS value\n(default: 0.8)')
        self.new_window.GroupA_Number_Edit.setToolTip('the number of sample in group A\n(default: 1)')
        self.new_window.GroupB_Number_Edit.setToolTip('the number of sample in group B\n(default: 1)')
        self.setWindowTitle(tool_name + ' ' + project_version)
        self.setFixedSize(self.width(), self.height())

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
        self.new_window.GroupA_FASTAQ_Path_Edit.setText(self.projectfile.FASTAQ_path_A)
        self.new_window.GroupB_FASTAQ_Path_Edit.setText(self.projectfile.FASTAQ_path_B)
        self.new_window.K_Value_Edit.setText(str(self.projectfile.K_value))
        self.new_window.CI_Value_Edit.setText(str(self.projectfile.Ci_value))
        self.new_window.CS_Value_Edit.setText(str(self.projectfile.Cs_value))
        self.new_window.KMC_Result_Path_Edit.setText(self.projectfile.KMC_path)
        self.new_window.Process_Number_Edit.setText(str(self.projectfile.Process_value))
        self.new_window.GM_Result_Path_Edit.setText(self.projectfile.GM_path)
        self.new_window.P_Value_Edit.setText(str(self.projectfile.P_value))
        self.new_window.ASS_n_Value_Edit.setText(str(self.projectfile.ASS_n_value))
        self.new_window.GF_Result_Path_Edit.setText(self.projectfile.GF_path)
        self.new_window.ASS_l_Value_Edit.setText(str(self.projectfile.ASS_l_value))
        self.new_window.GroupA_Number_Edit.setText(str(self.projectfile.GroupA_Number))
        self.new_window.GroupB_Number_Edit.setText(str(self.projectfile.GroupB_Number))

    def NewProject_Clicked(self):
        getdir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select the folder of a new project', './')
        if getdir == '': return
        self.project_dir = getdir
        if system_platform == 'Windows':
            self.project_dir = self.project_dir.replace('/', '\\')
        elif system_platform == 'Linux':
            self.project_dir = self.project_dir.replace('\\', '/')
        if os.listdir(self.project_dir) != []:
            QtWidgets.QMessageBox.critical(self, 'Error', 'Please select an empty folder!', QtWidgets.QMessageBox.Yes)
            return
        self.setWindowTitle(tool_name + ' ' + project_version + ' Project:' + self.project_dir)
        os.mkdir(os.path.join(self.project_dir, 'kmer_countings'))
        os.mkdir(os.path.join(self.project_dir, 'kmer_matrix'))
        os.mkdir(os.path.join(self.project_dir, 'kmer_features'))
        os.mkdir(os.path.join(self.project_dir, 'contig_result'))
        self.projectfile.CreateNewFile(self.project_dir)
        self.ReadConfiguration()

    def OpenProject_Clicked(self):
        getdir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select the floder of a project', './')
        if getdir == '': return
        self.project_dir = getdir
        if system_platform == 'Windows':
            self.project_dir = self.project_dir.replace('/', '\\')
        elif system_platform == 'Linux':
            self.project_dir = self.project_dir.replace('\\', '/')
        if os.path.exists(os.path.join(self.project_dir, 'ProjectList.list')) == False:
            QtWidgets.QMessageBox.critical(self, 'Error', 'Can not find a project file!', QtWidgets.QMessageBox.Yes)
            return
        self.setWindowTitle(tool_name + ' ' + project_version + ' Project:' + self.project_dir)
        self.projectfile.ReadFile(self.project_dir)
        self.ReadConfiguration()

    def SaveProject_Clicked(self):
        if self.project_dir == '':
            QtWidgets.QMessageBox.critical(self, 'Error', 'No project was opened!', QtWidgets.QMessageBox.Yes)
            return
        self.projectfile.WriteFile(self.project_dir)
        QtWidgets.QMessageBox.information(self, 'Success', 'Save this project successfully!', QtWidgets.QMessageBox.Yes)

    def Open_GroupA_FASTAQ_Button_Clicked(self):
        getdir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select the folder of A group sequencing files',
                                                            self.projectfile.FASTAQ_path_A)
        if getdir == '': return
        if system_platform == 'Windows':
            getdir = getdir.replace('/', '\\')
        elif system_platform == 'Linux':
            getdir = getdir.replace('\\', '/')
        self.projectfile.FASTAQ_path_A = getdir
        self.projectfile.GroupA_Number = len(os.listdir(getdir))
        self.new_window.GroupA_FASTAQ_Path_Edit.setText(self.projectfile.FASTAQ_path_A)
        self.new_window.GroupA_Number_Edit.setText(str(self.projectfile.GroupA_Number))

    def Open_GroupB_FASTAQ_Button_Clicked(self):
        getdir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select the folder of A group sequencing files',
                                                            self.projectfile.FASTAQ_path_B)
        if getdir == '': return
        if system_platform == 'Windows':
            getdir = getdir.replace('/', '\\')
        elif system_platform == 'Linux':
            getdir = getdir.replace('\\', '/')
        self.projectfile.FASTAQ_path_B = getdir
        self.projectfile.GroupB_Number = len(os.listdir(getdir))
        self.new_window.GroupB_FASTAQ_Path_Edit.setText(self.projectfile.FASTAQ_path_B)
        self.new_window.GroupB_Number_Edit.setText(str(self.projectfile.GroupB_Number))

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
        try: self.projectfile.ASS_n_value = float(self.new_window.ASS_n_Value_Edit.text())
        except: pass
        self.new_window.ASS_n_Value_Edit.setText(str(self.projectfile.ASS_n_value))

    def GroupA_Number_Edit_TextChange(self):
        try: self.projectfile.GroupA_Number = int(self.new_window.GroupA_Number_Edit.text())
        except: pass
        self.new_window.GroupA_Number_Edit.setText(str(self.projectfile.GroupA_Number))

    def GroupB_Number_Edit_TextChange(self):
        try: self.projectfile.GroupB_Number = int(self.new_window.GroupB_Number_Edit.text())
        except: pass
        self.new_window.GroupB_Number_Edit.setText(str(self.projectfile.GroupB_Number))

    def KMC_GO_Button_Clicked(self):
        if self.GO_flag == True:
            QtWidgets.QMessageBox.critical(self, 'Error', 'An assignment is running!', QtWidgets.QMessageBox.Yes)
            return
        self.GO_flag = True
        self.new_window.KMC_GO_Button.setEnabled(False)
        self.new_window.Open_GroupA_FASTAQ_Button.setEnabled(False)
        self.new_window.Open_GroupB_FASTAQ_Button.setEnabled(False)
        self.new_window.KMC_Result_Path_Button.setEnabled(False)
        self.new_window.K_Value_Edit.setEnabled(False)
        self.new_window.CI_Value_Edit.setEnabled(False)
        self.new_window.CS_Value_Edit.setEnabled(False)
        self.kmc_thread = kmc_read.KMC_Thread((self.projectfile.K_value, self.projectfile.Ci_value,
                                               self.projectfile.Cs_value, self.projectfile.FASTAQ_path_A,
                                               self.projectfile.FASTAQ_path_B, self.projectfile.KMC_path))
        self.kmc_thread.start()
        self.kmc_timer = QtCore.QTimer(self)
        self.kmc_timer.timeout.connect(self.KMC_Timer_Show)
        self.kmc_timer.start(100)

    def KMC_Timer_Show(self):
        if self.kmc_thread.status <= 0:
            self.GO_flag = False
            self.new_window.KMC_GO_Button.setEnabled(True)
            self.new_window.Open_GroupA_FASTAQ_Button.setEnabled(True)
            self.new_window.Open_GroupB_FASTAQ_Button.setEnabled(True)
            self.new_window.KMC_Result_Path_Button.setEnabled(True)
            self.new_window.K_Value_Edit.setEnabled(True)
            self.new_window.CI_Value_Edit.setEnabled(True)
            self.new_window.CS_Value_Edit.setEnabled(True)
            if self.kmc_thread.status == 0:
                self.projectfile.KMC_OK = True
                self.new_window.KMC_Project_Status_Label.setText('Status: Complete')
                self.new_window.KMC_Project_Status_Label.setStyleSheet('color:green')
            else:
                self.new_window.KMC_Project_Status_Label.setText('Status: ' + self.kmc_thread.loginfo)
                self.new_window.KMC_Project_Status_Label.setStyleSheet('color:red')
            del self.kmc_thread     # free up memery
            self.kmc_timer.stop()
        else:
            self.new_window.KMC_Project_Status_Label.setText('Status: ' + self.kmc_thread.loginfo)
            self.new_window.KMC_Project_Status_Label.setStyleSheet('color:blue')

    def GM_GO_Button_Clicked(self):
        if self.GO_flag == True:
            QtWidgets.QMessageBox.critical(self, 'Error', 'An assignment is running!', QtWidgets.QMessageBox.Yes)
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
        self.GO_flag = True
        self.new_window.KMC_Result_Path_Button.setEnabled(False)
        self.new_window.GM_Result_Path_Button.setEnabled(False)
        self.new_window.Process_Number_Edit.setEnabled(False)
        self.new_window.K_Value_Edit.setEnabled(False)
        self.new_window.CS_Value_Edit.setEnabled(False)
        self.new_window.GM_GO_Button.setEnabled(False)
        self.new_window.GroupA_Number_Edit.setEnabled(False)
        self.new_window.GroupB_Number_Edit.setEnabled(False)
        self.gm_thread = kmer_matrix.GM_Thread((self.projectfile.KMC_path, self.projectfile.GM_path,
                                               self.projectfile.Process_value, self.projectfile.K_value,
                                               self.projectfile.Cs_value, self.projectfile.GroupA_Number, self.projectfile.GroupB_Number))
        self.gm_thread.start()
        self.gm_timer = QtCore.QTimer(self)
        self.gm_timer.timeout.connect(self.GM_Timer_Show)
        self.gm_timer.start(1000)

    def GM_Timer_Show(self):
        if self.gm_thread.status <= 0:
            self.GO_flag = False
            self.new_window.KMC_Result_Path_Button.setEnabled(True)
            self.new_window.GM_Result_Path_Button.setEnabled(True)
            self.new_window.Process_Number_Edit.setEnabled(True)
            self.new_window.K_Value_Edit.setEnabled(True)
            self.new_window.CS_Value_Edit.setEnabled(True)
            self.new_window.GM_GO_Button.setEnabled(True)
            self.new_window.GroupA_Number_Edit.setEnabled(True)
            self.new_window.GroupB_Number_Edit.setEnabled(True)
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
            del self.gm_thread  # free up memery
            self.gm_timer.stop()
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
            QtWidgets.QMessageBox.critical(self, 'Error', 'An assignment is running!', QtWidgets.QMessageBox.Yes)
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
        self.GO_flag = True
        self.new_window.GF_GO_Button.setEnabled(False)
        self.new_window.KMC_Result_Path_Button.setEnabled(False)
        self.new_window.GM_Result_Path_Button.setEnabled(False)
        self.new_window.GF_Result_Path_Button.setEnabled(False)
        self.new_window.GroupA_Number_Edit.setEnabled(False)
        self.new_window.GroupB_Number_Edit.setEnabled(False)
        self.new_window.ASS_l_Value_Edit.setEnabled(False)
        self.new_window.P_Value_Edit.setEnabled(False)
        self.new_window.ASS_n_Value_Edit.setEnabled(False)
        self.gf_thread = kmer_features.GF_Thread((self.projectfile.GM_path, self.projectfile.GF_path,
                                                  self.projectfile.ASS_l_value, self.projectfile.P_value,
                                                  self.projectfile.ASS_n_value, self.projectfile.GroupA_Number,
                                                  self.projectfile.GroupB_Number))
        self.gf_thread.start()
        self.gf_timer = QtCore.QTimer(self)
        self.gf_timer.timeout.connect(self.GF_Timer_Show)
        self.gf_timer.start(1000)

    def GF_Timer_Show(self):
        if self.gf_thread.status <= 0:
            self.GO_flag = False
            self.new_window.GF_GO_Button.setEnabled(True)
            self.new_window.KMC_Result_Path_Button.setEnabled(True)
            self.new_window.GM_Result_Path_Button.setEnabled(True)
            self.new_window.GF_Result_Path_Button.setEnabled(True)
            self.new_window.GroupA_Number_Edit.setEnabled(True)
            self.new_window.GroupB_Number_Edit.setEnabled(True)
            self.new_window.ASS_l_Value_Edit.setEnabled(True)
            self.new_window.P_Value_Edit.setEnabled(True)
            self.new_window.ASS_n_Value_Edit.setEnabled(True)
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
            del self.gf_thread  # free up memery
            self.gf_timer.stop()
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
            QtWidgets.QMessageBox.critical(self, 'Error', 'An assignment is running!', QtWidgets.QMessageBox.Yes)
            return
        # initialization folders
        try:
            os.remove('temp')
        except:
            pass
        self.GO_flag = True
        self.new_window.KA_GO_Button.setEnabled(False)
        self.new_window.GF_Result_Path_Button.setEnabled(False)
        self.ka_thread = sequence_assembly.KA_Thread((self.projectfile.GF_path,
                                                      os.path.join(self.project_dir, 'contig_result')))
        self.ka_thread.start()
        self.ka_timer = QtCore.QTimer(self)
        self.ka_timer.timeout.connect(self.KA_Timer_Show)
        self.ka_timer.start(100)

    def KA_Timer_Show(self):
        if self.ka_thread.status <= 0:
            self.GO_flag = False
            self.new_window.KA_GO_Button.setEnabled(True)
            self.new_window.GF_Result_Path_Button.setEnabled(True)
            if self.ka_thread.status == 0:
                self.projectfile.KA_OK = True
                self.new_window.KA_Project_Status_Label.setText('Status: Complete')
                self.new_window.KA_Project_Status_Label.setStyleSheet('color:green')
            else:
                self.new_window.KA_Project_Status_Label.setText('Status: ' + self.ka_thread.loginfo)
                self.new_window.KA_Project_Status_Label.setStyleSheet('color:red')
            del self.ka_thread     # free up memery
            self.ka_timer.stop()
        else:
            self.new_window.KA_Project_Status_Label.setText('Status: ' + self.ka_thread.loginfo)
            self.new_window.KA_Project_Status_Label.setStyleSheet('color:blue')

if __name__ == "__main__":
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app = QtWidgets.QApplication(sys.argv)
    mywindow = myWindow()
    mywindow.show()
    sys.exit(app.exec_())
