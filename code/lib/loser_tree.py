# coding = utf-8
# author: QiChen
# version: v3.1
# modification date: 2020/3/20

from ctypes import *
import os, sys
import platform

class LoserTree(object):
    def __init__(self, kmerlist, KLEN, CLEN):
        self.KLEN = KLEN
        self.CLEN = CLEN
        system_platform = platform.system()
        if system_platform == 'Windows':
            self.lib = CDLL(os.path.join(os.getcwd(), r'.\lib\Loser_Tree_Lib.dll'))
        elif system_platform == 'Linux':
            file_path = os.path.dirname(os.path.realpath(sys.argv[0]))
            self.lib = CDLL(os.path.join(file_path, './lib/Loser_Tree_Lib.so'))
        self.lib.Loser_Tree_GetMin.argtypes = []
        self.lib.Loser_Tree_GetMin.restype = c_int
        self.kmer_list = ((c_char * (self.KLEN + 1)) * (self.CLEN + 1))()
        for i in range(len(kmerlist)):
            self.set_kmerlist(i, kmerlist[i])

    def set_kmerlist(self, lindex, kmer_fre):
        self.kmer_list[lindex].value = bytes(kmer_fre[:self.KLEN], 'utf-8')

    def setup_param(self):
        self.lib.Setup_parameters(self.KLEN, self.CLEN)

    def free_mem(self):
        self.lib.Free_memory()

    def adjust(self, tindex):
        self.lib.Loser_Tree_Adjust(tindex, pointer(self.kmer_list))

    def build(self):
        self.lib.Loser_Tree_Build(pointer(self.kmer_list))

    def getmin(self):
        return self.lib.Loser_Tree_GetMin()
