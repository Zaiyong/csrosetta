#!/usr/bin/env python2.7

# This function is used to get the noesyfile name in the peakfile


class Peaklist:
    def __init__(self,filelist,noesyfile_name):
        peaks=[]
        tmp_peaks=[]
        peak_cut_line=[]
        for i in range (0,len(noesyfile_name)):
            peak_cut_line.append(noesyfile_name[i][1])
        peak_cut_line.append(len(filelist))
        for i in range (0,len(noesyfile_name)):
            for j in range(peak_cut_line[i],peak_cut_line[i+1]):
                if filelist[j].split()[0][0]<='9' and filelist[j].split()[0][0]>='0':
                    tmp_peaks.append([filelist[j].split()[0],filelist[j].split()[1:noesyfile_name[i][2]+1],filelist[j].split()[noesyfile_name[i][2]+3]])
            peaks.append(tmp_peaks)
            tmp_peaks=[]
    
        self.data=peaks
        

