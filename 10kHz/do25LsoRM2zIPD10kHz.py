# -*- coding: utf-8 -*-
"""
do25LsoRM2zIPD10kHz.py

Calls a model LSO neuron type "2z" (brisk membrane, lson from
lson25RM2z.py) with input from spikeTime files of model auditory nerve fibers 
(Zilany et al., 2014) for 25 stimulus repetitions of transposed tones at 10 kHz,
with modulation rates given in variable modFreqs_Hz.
Input files are drawn from, and output files are written to, subfolders:
    stimReps\rep## where ## ranges from 01 to 25.

Created on Fri 31 May 2019

@author: Andrew Brughera
"""
import os
from lson25RM2z import lson

NStimReps = 25
modFreqs_Hz = [32, 64, 128, 256, 512, 800]
EPhsStep_deg = 15 # EPhs: envelope phase (of the modulation period)

# Commenting the next line to Assume the Initial Directory  
#os.chdir('C:\\Users\\Andrew Brughera\\lso\\IPD\\10kHz')
# will allow for an error message from gbCoSpkFiles.append('ANSpTsCaCF10kHz...
# if we don't start in the 10kHz directory.
os.chdir('stimReps')
 
for stimRep in range(NStimReps):
#for stimRep in range(2,NStimReps): # (modify to continue a partial run)
    stimRepPlus1Str = str(stimRep+1).zfill(len(str(NStimReps)))
    os.chdir('rep' + stimRepPlus1Str) # 'rep01'..'rep25' for Matlab
    
    for modF in modFreqs_Hz:
        # As in Wang & Colburn (2012), the LSO model has simplified inputs.
        # Spike times from model auditory nerve fibers represent inputs
        # driven by the cochlear nucleus: contralateral inhibitory inputs
        # driven by globular bushy cells via the MNTB, and ipsilateral
        # excitatory inputs from spherical bushy cells.
        gbCoSpkFiles = [] # Lists of string filenames are named "---Files". 
        sbIpSpkFiles = []
        
        # EPhs: envelope phase (of the modulation period)
        # EIPD: envelope interaural-phase-difference
        # 24 instances of 25 EIPDs, -180:15:180.
        # For EIPD +180 (Contralateral-leading):
        # gbCo EPhs 0:-15:-345; sbIp EPhs -180:-15:-525
        for EPhsIdx in range(int(540/EPhsStep_deg)):
            gbCoSpkFiles.append('ANSpTsCaCF10kHz1TT'+ str(modF).zfill(3) + 'HzEPhsN'+ str(EPhsIdx*EPhsStep_deg).zfill(3) + 'SPL75dB' + stimRepPlus1Str + 'MdSp3.txt')
            sbIpSpkFiles.append('ANSpTsCaCF10kHz2TT'+ str(modF).zfill(3) + 'HzEPhsN'+ str(EPhsIdx*EPhsStep_deg).zfill(3) + 'SPL75dB' + stimRepPlus1Str + 'MdSp3.txt')
        
        # 13 ipsi-leading EIPDs 0, -15, -30, -45 ... -180 degrees:
        for EIPDIdx in range(0,int((180/EPhsStep_deg)+1)):
            #lsonOutputTuples = []
            # 24 bilateral pairs of inputs: ipsi EPhs 0:-15:-345
            for ipEPhsIdx in range(int(360/EPhsStep_deg)):
                # Note, in LSO model contra is 1st, Ipsi is 2nd
                inputSpkFileTuple = (gbCoSpkFiles[ipEPhsIdx+EIPDIdx], sbIpSpkFiles[ipEPhsIdx])
                lsonOutputTuple = lson(inputSpkFileTuple)
                del lsonOutputTuple # reduce memory used
                #lsonOutputTuples.append(lsonOutputTuple)
            # end of lson loop
        # end of ITD shift loop
        # 12 contra-leading EIPDs +15, +30, +45 ... +180 degrees:
        for EIPDIdx in range(1,int((180/EPhsStep_deg)+1)):
            #lsonOutputTuples = []
            # 24 bilateral pairs of inputs: contra EPhs 0:-15:-345
            for coEPhsIdx in range(int(360/EPhsStep_deg)):
                # Note, in LSO model contra is 1st, Ipsi is 2nd
                inputSpkFileTuple = (gbCoSpkFiles[coEPhsIdx], sbIpSpkFiles[coEPhsIdx+EIPDIdx])
                lsonOutputTuple = lson(inputSpkFileTuple)
                del lsonOutputTuple # reduce memory used
                #lsonOutputTuples.append(lsonOutputTuple)
            # end of lson loop
        # end of ITD shift loop
    # end of modF loop
    #os.chdir('C:\\Users\\Andrew Brughera\\lso\\IPD\\10kHz\\stimReps')
    os.chdir('..')
# end of stimRep loop
os.chdir('..')
os.getcwd() # should return 'C:\\Users\\Andrew Brughera\\lso\\IPD\\10kHz'
