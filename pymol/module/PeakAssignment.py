#!/usr/bin/env python2.7

# This function is used to get the noesyfile name in the peakfile


class PeakAssignment:
    def __init__(self,input):
        self.data=input
    filelist=silf.data
    noesyfile_name=NOESY_Data_collection(filelist).data
    peakassign=[]
    tmp_peakassign=[]
    peak_cut_line=[]
    for i in range (0,len(noesyfile_name)):
        peak_cut_line.append(noesyfile_name[i][1])
    peak_cut_line.append(len(filelist))
    for i in range (0,len(noesyfile_name)):
        for j in range(peak_cut_line[i],peak_cut_line[i+1]):
            if '#W' in filelist[j].split() and ( '#eliminated:' not in filelist[j].split()):
                if filelist[j].split()[0][0]<='9' and filelist[j].split()[0][0]>='0':
                    tmp_index=filelist[j].split()[0]
                    tmp_peak=filelist[j].split()[noesyfile_name[i][2]+3]
                    e_index=filelist[j].split().index('e')
                    W_index=filelist[j].split().index('#W')
                    tmp_peakassign.append([tmp_index,tmp_peak,filelist[j].split()[e_index+2:W_index]])
                else:
                    W_index=filelist[j].split().index('#W')
                    tmp_peakassign.append([tmp_index,tmp_peak,filelist[j].split()[0:W_index]])

        
        peakassign.append(tmp_peakassign)
        tmp_peakassign=[]
    
    self.data=peakassign
        

