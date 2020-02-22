# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 16:18:14 2018

@author: Andrew Brughera
"""


#files = [] # list of tuple with old filename and new filename
#for file in os.listdir(path):
#    if fname.startswith('AN'):
#        if file.endswith(".txt"):
#            if file.find("Carrier") > 0:
#                counter = counter + 1
#                newFileName = file.replace("_ready", "_busy"))
#                file = newFileName
#                files.append((file, newFileName))
            
            
            
from os import rename, listdir

#fnames = listdir('.')
#for fname in fnames:
#    if fname.startswith('Lso'):
#        if fname.endswith('.txt'):
#            newFileName = fname.replace('Ge10Te0p2Gi60Ti1', 'Te0p2Ge10Ti1Gi60')
#            rename(fname, newFileName)
            
#fnames = listdir('.')
#for fname in fnames:
#    if fname.startswith('Lso'):
#        if fname.endswith('.txt'):
#            fname1 = fname.replace('Spec3.', 'Spec3Ge10Te0p2Gi60Ti1.')
#            fname2 = fname1.replace('AM64Hz', 'AM064Hz')
#            fname3 = fname2.replace('Phs000', 'PhsN000')
#            fname4 = fname3.replace('CF04kHzCa04kHz', 'CaCF04kHz')
#            fname5 = fname4.replace('CF06kHzCa06kHz', 'CaCF06kHz')
#            fname6 = fname5.replace('CF10kHzCa10kHz', 'CaCF10kHz')
#            fname7 = fname6.replace('SpkTmSec', 'SpTs')
#            newFileName = fname7.replace('EnvPhs', 'EPhs')
#            rename(fname, newFileName) 
            
            
#fnames = listdir('.')
#for fname in fnames:
#    if fname.startswith('Lso'):
#        if fname.endswith('.txt'):
#            fname1 = fname.replace('IPD0C', 'IPD000C')
#            fname2 = fname1.replace('IPD15C', 'IPD015C')
#            fname3 = fname2.replace('IPD30C', 'IPD030C')
#            fname4 = fname3.replace('IPD45C', 'IPD045C')
#            fname5 = fname4.replace('IPD60C', 'IPD060C')
#            fname6 = fname5.replace('IPD75C', 'IPD075C')
#            fname7 = fname6.replace('IPD90C', 'IPD090C')
#            fname8 = fname7.replace('IPDN15C', 'IPDN015C')
#            fname9 = fname8.replace('IPDN30C', 'IPDN030C')
#            fname10 = fname9.replace('IPDN45C', 'IPDN045C')
#            fname11 = fname10.replace('IPDN60C', 'IPDN060C')
#            fname12 = fname11.replace('IPDN75C', 'IPDN075C')
#            newFileName = fname12.replace('IPDN90C', 'IPDN090C')
#            rename(fname, newFileName)
            
fnames = listdir('.')
for fname in fnames:
    if fname.startswith('Lso'):
        if fname.endswith('.txt'):
            fname1 = fname.replace('IPD0', 'IPDP0')
            fname2 = fname1.replace('IPD1', 'IPDP1')
            fname3 = fname2.replace('IPD2', 'IPDP2')
            fname4 = fname3.replace('IPD3', 'IPDP3')
            fname5 = fname4.replace('IPD4', 'IPDP4')
            fname6 = fname5.replace('IPD5', 'IPDP5')
            fname7 = fname6.replace('IPD6', 'IPDP6')
            fname8 = fname7.replace('IPD7', 'IPDP7')
            fname9 = fname8.replace('IPD8', 'IPDP8')
            newFileName = fname9.replace('IPD9', 'IPDP9')
            rename(fname, newFileName)
