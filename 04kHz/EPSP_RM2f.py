# -*- coding: utf-8 -*-
"""
Downloaded from Brian2 website on Sat Jan 13 10:02:04 2018

Cochlear neuron model of Rothman & Manis
----------------------------------------
Rothman JS, Manis PB (2003) The roles potassium currents play in
regulating the electrical activity of ventral cochlear nucleus neurons.
J Neurophysiol 89:3097-113.

All model types differ only by the maximal conductances.

Adapted from their Neuron implementation by Romain Brette
"""
#from brian2 import *
from brian2 import mV, pA, pF, nS, ms, second, Hz
from brian2 import SpikeGeneratorGroup, NeuronGroup, Synapses
from brian2 import SpikeMonitor, StateMonitor, run
from brian2 import defaultclock
from math import exp
import numpy as np
import matplotlib.pyplot as plt
    
defaultclock.dt=0.02*ms # for better precision

'''
Simulation parameters: choose current amplitude and neuron type
(from type1c, type1t, type12, type 21, type2, type2o)
'''
neuron_type = 'type2g2x'
nLsons = 1 # number of LSO neurons
temp_degC=37.
Vrest = -63.6*mV # observed resting potential for model 2f

#    nAnfsPerInputFile = 40
#
#    nGbcsCoPerLson = 8  # Gjoni et al. 2018
#    nSbcsIpPerLson = 40 # Gjoni et al. 2018
#
#    anCoSpkFile = lsoInputSpkFileTuple[0]
#    anIpSpkFile = lsoInputSpkFileTuple[1]
#
#    gbCoIdxSpktimeArray = np.loadtxt(anCoSpkFile)
#    gbCoCellIndices = gbCoIdxSpktimeArray[:, 0].astype(int)
#    # For now, spiketimes from Zilany AN in SECONDS, so * 1000*ms
#    gbCoSpkTimes = gbCoIdxSpktimeArray[:, 1] * 1000. * ms
#    gbCoSpkGenGrp = SpikeGeneratorGroup(nAnfsPerInputFile, gbCoCellIndices, gbCoSpkTimes)
#    
#    sbIpIdxSpktimeArray = np.loadtxt(anIpSpkFile)
#    sbIpCellIndices = sbIpIdxSpktimeArray[:, 0].astype(int)
#    # For now, spiketimes from Zilany AN in SECONDS, so * 1000*ms
#    sbIpSpkTimes = sbIpIdxSpktimeArray[:, 1] * 1000. * ms
#    sbIpSpkGenGrp = SpikeGeneratorGroup(nAnfsPerInputFile, sbIpCellIndices, sbIpSpkTimes)

nAnfsPerInputFile = 40

nGbcsCoPerLson = 8  # Gjoni et al. 2018
nSbcsIpPerLson = 40 # Gjoni et al. 2018
#nGbcsCoPerLson = 1  # showing unitary IPSP
#nSbcsIpPerLson = 1 # showing unitary EPSP

#gbCoCellIndices = np.array([0])
#gbCoSpkTimes = np.array([10.0]) * ms

sbIpCellIndices = np.array([0])
sbIpSpkTimes = np.array([5.0]) * ms

#gbCoSpkGenGrp = SpikeGeneratorGroup(nAnfsPerInputFile, gbCoCellIndices, gbCoSpkTimes)
sbIpSpkGenGrp = SpikeGeneratorGroup(nAnfsPerInputFile, sbIpCellIndices, sbIpSpkTimes)

C = 12*pF
EH = -43*mV
EK = -70*mV  # -77*mV in mod file
EL = -65*mV
ENa = 55*mV # 55*mV in RM2003; 50*mv by Brette
nf = 0.85  # proportion of n vs p kinetics
zss = 0.5  # steady state inactivation of glt
# default temp_degC = 37., human body temperature in degree celcius
# q10 for ion-channel time-constants (RM2003, p.3106):
q10 = 3. ** ((temp_degC - 22.) / 10.)
# q10 for ion-channel gbar parameters (RM2003, p.3106):
q10gbar = 2. ** ((temp_degC - 22.) / 10.)
# hcno current (octopus cell)
frac = 0.0
qt = 4.5 ** ((temp_degC - 33.) / 10.)
# Synaptic parameters:
Es_e = 0.*mV
tausE = 0.2*ms
Es_i = -90*mV
tausI = 0.4*ms
'''Synaptic weights are unitless according to Brian2.
The effective unit is siemens, so they can work in amp, volt, siemens eqns.
We multiply synaptic weight w_e by unit siemens when adding it to g_e.
We use a local variable w_e for synaptic weight rather than the standard w:''' 
w_elson = 10e-9
w_ilson = 100e-9 # 6e-9 @ 200Hz; 12e-9 @ 600 Hz
'''Here's why:
The synapses sbc3SynE.w synaptic weight references the Same State Variable as 
as the neuron group sbc3Grp.w (klt activation w).
Update sbc3Grp.w, and you set sbc3SynE.w to same value, and vice versa.
This way klt activation and synaptic weight are identical and quite ridiculous.
So use a local variable other than w for the synaptic weight!'''

# Maximal conductances of different cell types in nS
maximal_conductances = dict(
type1c=(1000, 150, 0, 0, 0.5, 0, 2),
type1t=(1000, 80, 0, 65, 0.5, 0, 2),
type12=(1000, 150, 20, 0, 2, 0, 2),
type21=(1000, 150, 35, 0, 3.5, 0, 2),
type2=(1000, 150, 200, 0, 20, 0, 2),
type2g2x=(2000, 300, 400, 0, 40, 0, 2),
type2g1p5x=(1000, 150, 300, 0, 30, 0, 2),
type2g1p2x=(1200, 180, 240, 0, 24, 0, 2),
type2g0p5x=(1000, 150, 100, 0, 10, 0, 2),
type2o=(1000, 150, 600, 0, 0, 40, 2) # octopus cell
)
gnabar, gkhtbar, gkltbar, gkabar, ghbar, gbarno, gl = [x * nS for x in maximal_conductances[neuron_type]]

# Classical Na channel
eqs_na = """
ina = gnabar*m**3*h*(ENa-v) : amp
dm/dt=q10*(minf-m)/mtau : 1
dh/dt=q10*(hinf-h)/htau : 1
minf = 1./(1+exp(-(vu + 38.) / 7.)) : 1
hinf = 1./(1+exp((vu + 65.) / 6.)) : 1
mtau =  ((10. / (5*exp((vu+60.) / 18.) + 36.*exp(-(vu+60.) / 25.))) + 0.04)*ms : second
htau =  ((100. / (7*exp((vu+60.) / 11.) + 10.*exp(-(vu+60.) / 25.))) + 0.6)*ms : second
"""

# KHT channel (delayed-rectifier K+)
eqs_kht = """
ikht = gkhtbar*(nf*n**2 + (1-nf)*p)*(EK-v) : amp
dn/dt=q10*(ninf-n)/ntau : 1
dp/dt=q10*(pinf-p)/ptau : 1
ninf = (1 + exp(-(vu + 15) / 5.))**-0.5 : 1
pinf =  1. / (1 + exp(-(vu + 23) / 6.)) : 1
ntau = ((100. / (11*exp((vu+60) / 24.) + 21*exp(-(vu+60) / 23.))) + 0.7)*ms : second
ptau = ((100. / (4*exp((vu+60) / 32.) + 5*exp(-(vu+60) / 22.))) + 5)*ms : second
"""

# Ih channel (subthreshold adaptive, non-inactivating)
eqs_ih = """
ih = ghbar*r*(EH-v) : amp
dr/dt=q10*(rinf-r)/rtau : 1
rinf = 1. / (1+exp((vu + 76.) / 7.)) : 1
rtau = ((100000. / (237.*exp((vu+60.) / 12.) + 17.*exp(-(vu+60.) / 14.))) + 25.)*ms : second
"""

# KLT channel (low threshold K+)
eqs_klt = """
iklt = gkltbar*w**4*z*(EK-v) : amp
dw/dt=q10*(winf-w)/wtau : 1
dz/dt=q10*(zinf-z)/ztau : 1
winf = (1. / (1 + exp(-(vu + 48.) / 6.)))**0.25 : 1
zinf = zss + ((1.-zss) / (1 + exp((vu + 71.) / 10.))) : 1
wtau = ((100. / (6.*exp((vu+60.) / 6.) + 16.*exp(-(vu+60.) / 45.))) + 1.5)*ms : second
ztau = ((1000. / (exp((vu+60.) / 20.) + exp(-(vu+60.) / 8.))) + 50)*ms : second
"""

# Ka channel (transient K+)
eqs_ka = """
ika = gkabar*a**4*b*c*(EK-v): amp
da/dt=q10*(ainf-a)/atau : 1
db/dt=q10*(binf-b)/btau : 1
dc/dt=q10*(cinf-c)/ctau : 1
ainf = (1. / (1 + exp(-(vu + 31) / 6.)))**0.25 : 1
binf = 1. / (1 + exp((vu + 66) / 7.))**0.5 : 1
cinf = 1. / (1 + exp((vu + 66) / 7.))**0.5 : 1
atau =  ((100. / (7*exp((vu+60) / 14.) + 29*exp(-(vu+60) / 24.))) + 0.1)*ms : second
btau =  ((1000. / (14*exp((vu+60) / 27.) + 29*exp(-(vu+60) / 24.))) + 1)*ms : second
ctau = ((90. / (1 + exp((-66-vu) / 17.))) + 10)*ms : second
"""

# Leak
eqs_leak = """
ileak = gl*(EL-v) : amp
"""

# h current for octopus cells
eqs_hcno = """
ihcno = gbarno*(h1*frac + h2*(1-frac))*(EH-v) : amp
dh1/dt=(hinfno-h1)/tau1 : 1
dh2/dt=(hinfno-h2)/tau2 : 1
hinfno = 1./(1+exp((vu+66.)/7.)) : 1
tau1 = bet1/(qt*0.008*(1+alp1))*ms : second
tau2 = bet2/(qt*0.0029*(1+alp2))*ms : second
alp1 = exp(1e-3*3*(vu+50)*9.648e4/(8.315*(273.16+temp_degC))) : 1
bet1 = exp(1e-3*3*0.3*(vu+50)*9.648e4/(8.315*(273.16+temp_degC))) : 1 
alp2 = exp(1e-3*3*(vu+84)*9.648e4/(8.315*(273.16+temp_degC))) : 1
bet2 = exp(1e-3*3*0.6*(vu+84)*9.648e4/(8.315*(273.16+temp_degC))) : 1
"""

eqs = """
dv/dt = (ileak + ina + ikht + iklt + ika + ih + ihcno + I + Is_e + Is_i)/C : volt
Is_e = gs_e * (Es_e - v) : amp
gs_e : siemens
Is_i = gs_i * (Es_i - v) : amp
gs_i : siemens
vu = v/mV : 1  # unitless v
I : amp
"""
#Added Is_i to RM2003

eqs += eqs_leak + eqs_ka + eqs_na + eqs_ih + eqs_klt + eqs_kht + eqs_hcno


lsonGrp = NeuronGroup(nLsons, eqs, method='exponential_euler', 
                    threshold='v > -30*mV', refractory='v > -45*mV')
#gbcGrp.I = 2500.0*pA
lsonGrp.I = 0.0*pA
# Initialize model near v_rest with no inputs
lsonGrp.v = Vrest
vu = lsonGrp.v/mV # unitless v
lsonGrp.m = 1./(1+exp(-(vu + 38.) / 7.))
lsonGrp.h = 1./(1+exp((vu + 65.) / 6.))
lsonGrp.n = (1 + exp(-(vu + 15) / 5.))**-0.5
lsonGrp.p = 1. / (1 + exp(-(vu + 23) / 6.))
lsonGrp.r = 1. / (1+exp((vu + 76.) / 7.))
lsonGrp.w = (1. / (1 + exp(-(vu + 48.) / 6.)))**0.25
lsonGrp.z = zss + ((1.-zss) / (1 + exp((vu + 71.) / 10.)))
lsonGrp.a = (1. / (1 + exp(-(vu + 31) / 6.)))**0.25
lsonGrp.b = 1. / (1 + exp((vu + 66) / 7.))**0.5
lsonGrp.c = 1. / (1 + exp((vu + 66) / 7.))**0.5
lsonGrp.h1 = 1./(1+exp((vu+66.)/7.))
lsonGrp.h2 = 1./(1+exp((vu+66.)/7.))


#lsonSynI = Synapses(gbCoSpkGenGrp, lsonGrp, 
#                 model='''dg_i/dt = -g_i/tausI : siemens (clock-driven)
#                          gs_i_post = g_i : siemens (summed)''',
#                 on_pre='g_i += w_ilson*siemens',
#                 method = 'exact')
#lsonSynI.connect(i=np.arange(nGbcsCoPerLson), j=0)

lsonSynE = Synapses(sbIpSpkGenGrp, lsonGrp, 
                 model='''dg_e/dt = -g_e/tausE : siemens (clock-driven)
                          gs_e_post = g_e : siemens (summed)''',
                 on_pre='g_e += w_elson*siemens',
                 method = 'exact')
lsonSynE.connect(i=np.arange(nSbcsIpPerLson), j=0)

lsonSpks = SpikeMonitor(lsonGrp)
lsonState = StateMonitor(lsonGrp, ['gs_e','Is_e','gs_i','Is_i','v'], record=True)

run(50*ms, report='text')

plt.plot(lsonState.t / ms, lsonState[0].gs_e / nS)
plt.xlabel('t (ms)')
plt.ylabel('gs_e (nS)')
plt.show()

plt.plot(lsonState.t / ms, lsonState[0].Is_e / pA)
plt.xlabel('t (ms)')
plt.ylabel('Is_e (pA)')
plt.show()

plt.plot(lsonState.t / ms, lsonState[0].gs_i / nS)
plt.xlabel('t (ms)')
plt.ylabel('gs_i (nS)')
plt.show()

plt.plot(lsonState.t / ms, lsonState[0].Is_i / pA)
plt.xlabel('t (ms)')
plt.ylabel('Is_i (pA)')
plt.show()

plt.plot(lsonState.t / ms, lsonState[0].v / mV)
plt.xlabel('t (ms)')
plt.ylabel('v (mV)')
plt.show()

# Synaptic parameters in output filename
if (abs(round(tausE/ms)-(tausE/ms)) < 0.1):
    Te = str(round(tausE/ms))
else:
    Te = str(round(10*tausE/ms)/10.)
Te = Te.replace('.','p')

if (abs(round(w_elson/1e-9)-(w_elson/1e-9)) < 0.1):
    We = str(round(w_elson/1e-9))
else:
    We = str(round(10*w_elson/1e-9)/10)
We = We.replace('.','p')

VmFileStr = 'EPSPmV_RM2f_' + 'tauE'+Te+'ms_' 'gE'+We+'nS' + '.txt'
file9 = open(VmFileStr,'w')
for index in range(len(lsonState.t)):
    file9.write(str(lsonState.t[index] / ms) + " " + str(lsonState[0].v[index] / mV) + '\n')
file9.close()
