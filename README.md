# LSO_Transposed_Tones
 (Anaconda Python 3.7, Brian 2 Neural Simulator)HodgkinHuxleyType models: EI neurons in populations

 Lateral Superior Olive (LSO) EI neurons:
 40 ipsilateral excitatory (E) inputs
 8 contralateral inhibitory (I) inputs
  
 Each input is driven by a model 
 auditory nerve fiber (ANF) (Zilany et al. 2014)
 of medium spontaneous rate,
 and high characteristic-frequency.
 Complete ANF Matlab code downloaded from
 https://www.urmc.rochester.edu/labs/carney/publications-code/auditory-models.aspx

 24 model LSO neurons per population
 6 populations by membrane properties:
 1c - Slow - Rothman & Manis 2003 (RM03);
 1a - Moderate;
 1b - HL-Quick;
 2z - Brisk;
 2  - Fast - Wang & Colburn 2012, RM03;
 2f - Very Fast.
 Brisk, fast & very fast include
 low-threshold potassium (KLT) currents.

 Transposed tone stimuli:
 Carrier frequencies: 4, 6, 10 kHz
 Modulation rates: 32, 64, 128, 256, 512, 800 Hz
 Envelope is half-wave rectified sinusoid,
 with no frequency components above 2 kHz.
 25 stimulus repetitions per condition:
 do25*.py call lson25*.py.

 In related filenames, IPD refers to envelope interaural phase difference, the IPD with respect to the amplitude modulation cycle. Phase in the carrier wareform was randomized in all stimuli at each envelope phase.
