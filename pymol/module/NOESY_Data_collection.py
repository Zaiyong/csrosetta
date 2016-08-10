#!/usr/bin/env python2.7

# This function is used to get the noesyfile name in the peakfile


class NOESY_Data_collection:
    def __init__(self,input):
        self.data=input
    filelist=input
    noesyfile_name=[]
    for i in range (0,len(filelist)):
        if filelist[i].split()[0]=='#FILENAME':
            tmp_ilename=filelist[i].split()[1]
            tmp_line=i
        if filelist[i].split()[0]=='#TOLERANCE':
            noesyfile_name.append([tmp_ilename,tmp_line,len(filelist[i].split())-1])
    def result(self):
        ressult.data=noesyfile_name
#return noesyfile_name
        

