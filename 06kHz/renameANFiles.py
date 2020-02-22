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

fnames = listdir('.')
for fname in fnames:
    if fname.startswith('AN'):
        if fname.endswith('.txt'):
            fname1 = fname.replace('TT32Hz', 'TT032Hz')
            fname2 = fname1.replace('Phs000', 'PhsN000')
            fname3 = fname2.replace('CF4kHzCa4kHz', 'CaCF04kHz')
            fname4 = fname3.replace('CF6kHzCa6kHz', 'CaCF06kHz')
            fname5 = fname4.replace('CF10kHzCa10kHz', 'CaCF10kHz')
            newFileName = fname5.replace('TT64Hz', 'TT064Hz')
            rename(fname, newFileName)