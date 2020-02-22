# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 18:38:23 2019

@author: Andrew Brughera
"""

import os
from os import rename, listdir

NStimReps = 25

os.chdir('stimReps\\OffCF')
 
for stimRep in range(NStimReps):
    stimRepPlus1Str = str(stimRep+1).zfill(len(str(NStimReps)))
    os.chdir('rep' + stimRepPlus1Str) # 'rep01'..'rep25' for Matlab
    fnames = listdir('.')
    for fname in fnames:
        if fname.startswith('Lso'):
            if fname.endswith('.txt'):
                newFileName = fname.replace('SpTms', 'SpT')
                newFileName = newFileName.replace('kCa', 'Ca')
                newFileName = newFileName.replace('Co', 'C')
                newFileName = newFileName.replace('EPhs', 'EP')
                rename(fname, newFileName)
    os.chdir('..')
                
os.chdir('..\\..')
