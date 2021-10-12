# -*- coding: utf-8 -*-
"""
doLsoRM2f_IID_4kHz.py

Calls a model LSO neuron type 2f (resonant, very fast membrane, 
lson from lsonRM2f_IID.py) with input from spikeTime files of model auditory 
nerve fibers (Zilany et al., 2014) for stimulus pure tones at 4 kHz.

Created on Wed 5 May 2021

@author: Andrew Brughera
"""
#import os
from lsonRM2f_IID import lson

dBstep = 2;
for SPLdB in range(0,61,dBstep):
    # inputSpkFileTuple[0] (Ear1) is contralateral; inputSpkFileTuple[1] (Ear2) is ipsilateral
    inputSpkFileTuple = ('ANSpTsToneCaCF04kHzEar1dBSPL' + str(SPLdB).zfill(2) + 'MdSp3.txt', 'ANSpTsToneCaCF04kHzEar2dBSPL30MdSp3.txt')
    lsonOutputTuple = lson(inputSpkFileTuple)
    
for SPLdB in range(10,71,dBstep):
    inputSpkFileTuple = ('ANSpTsToneCaCF04kHzEar1dBSPL' + str(SPLdB).zfill(2) + 'MdSp3.txt', 'ANSpTsToneCaCF04kHzEar2dBSPL40MdSp3.txt')
    lsonOutputTuple = lson(inputSpkFileTuple)

for SPLdB in range(20,81,dBstep):
    inputSpkFileTuple = ('ANSpTsToneCaCF04kHzEar1dBSPL' + str(SPLdB).zfill(2) + 'MdSp3.txt', 'ANSpTsToneCaCF04kHzEar2dBSPL50MdSp3.txt')
    lsonOutputTuple = lson(inputSpkFileTuple)
    
