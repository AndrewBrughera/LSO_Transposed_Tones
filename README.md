# LSO_Transposed_Tones
 HodgkinHuxleyType models: EI neurons in populations

 24 model neurons per population
 4 populations by speed:
 1c - Slow - Rothman & Manis 2003 (RM03)
 2z - Moderate
 2  - Fast - Wang & Colburn 2012, RM03
 2f - Very Fast
 Moderate, fast & very fast include
 low-threshold potassium (KLT) currents.

 Transposed tone stimuli:
 Carrier frequencies: 4, 6, 10 kHz
 Modulation rates: 32, 64, 128, 256, 512, 800 Hz
 Envelope is half-wave rectified sinusoid,
 with no frequency components above 2 kHz.
 25 stimulus repetitions per condition:
 do25*.py call lson25*.py.
