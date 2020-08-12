# -*- coding: utf-8 -*-
"""
do25CF4p28LsoRM2zIPD4kHz.py

Calls a model LSO neuron type 2z (brisk membrane)
with input from spikeTime files of model auditory nerve fibers 
(Zilany et al., 2014) for 25 stimulus repetitions of transposed tones at 4 kHz,
with modulation rates given in variable modFreqs_Hz.
Input files are drawn from, and output files are written to, subfolders:
    stimReps\rep## where ## ranges from 01 to 25.

Created on Sat 14 Sep 2019

@author: Andrew Brughera
"""
import os
from lson25CF4pXXRM2z import lson

NStimReps = 25
modFreqs_Hz = [512, 800]
EPhsStep_deg = 15 # EPhs: envelope phase (of the modulation period)

# Commenting the next line to Assume the Initial Directory  
#os.chdir('C:\\Users\\Andrew Brughera\\lso\\IPD\\4kHz')
# will allow for an error message from gbCoSpkFiles.append('ANSpTsCaCF04kHz...
# if we don't start in the 4kHz directory.
os.chdir('stimReps\\OffCF')
 
for stimRep in range(NStimReps):
#for stimRep in range(15,NStimReps):
    stimRepPlus1Str = str(stimRep+1).zfill(len(str(NStimReps)))
    os.chdir('rep' + stimRepPlus1Str) # 'rep01'..'rep25' for Matlab
    
    for modF in modFreqs_Hz:
        gbCoSpkFiles = []
        sbIpSpkFiles = []
        
        # EPhs: envelope phase (of the modulation period)
        # EIPD: envelope interaural-phase-difference
        # 24 instances of 25 EIPDs, -180:15:180.
        # For EIPD +180 (Contralateral-leading):
        # gbCo EPhs 0:-15:-345; sbIp EPhs -180:-15:-525
        # ANSpTsCF04p28kHzCa04kHz1TT512HzEPhsN000SPL75dB18MdSp3
        for EPhsIdx in range(int(540/EPhsStep_deg)):
            gbCoSpkFiles.append('ANSpTsCF04p28kHzCa04kHz1TT'+ str(modF).zfill(3) + 'HzEPhsN'+ str(EPhsIdx*EPhsStep_deg).zfill(3) + 'SPL75dB' + stimRepPlus1Str + 'MdSp3.txt')
            sbIpSpkFiles.append('ANSpTsCF04p28kHzCa04kHz2TT'+ str(modF).zfill(3) + 'HzEPhsN'+ str(EPhsIdx*EPhsStep_deg).zfill(3) + 'SPL75dB' + stimRepPlus1Str + 'MdSp3.txt')
        
        # 13 ipsi-leading IPDs 0, -15, -30, -45 ... -180 degrees:
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
        # 12 contra-leading ITDs +15, +30, +45 ... +180 degrees:
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
    #os.chdir('C:\\Users\\Andrew Brughera\\lso\\IPD\\4kHz\\stimReps')
    os.chdir('..')
# end of stimRep loop
os.chdir('..')
os.getcwd() # should return 'C:\\Users\\Andrew Brughera\\lso\\IPD\\4kHz'
