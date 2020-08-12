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
from brian2 import NeuronGroup
from brian2 import StateMonitor, run
from brian2 import defaultclock
from math import exp
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.ticker import ScalarFormatter
#for axis in [ax.xaxis, ax.yaxis]:
#    axis.set_major_formatter(ScalarFormatter())
    
defaultclock.dt=0.02*ms # for better precision

'''
Simulation parameters: choose current amplitude and neuron type
(from type1c, type1t, type12, type 21, type2, type2o)
'''
neuron_type = 'LsoNp1a'
nLsons = 1 # number of LSO neurons
temp_degC=22.
Vrest = -64.1*mV # observed resting potential for model LsoNp1a

C = 21.7*pF # 21.7pF: non-principal LSO neuron (Barnes-Davies et al. 2004)
Eh = -43*mV
EK = -70*mV  # -77*mV in mod file
El = -65*mV
ENa = 50*mV
nf = 0.85  # proportion of n vs p kinetics
zss = 0.5  # steady state inactivation of glt
q10 = 3. ** ((temp_degC - 22) / 10.)
# hcno current (octopus cell)
frac = 0.0
qt = 4.5 ** ((temp_degC - 33.) / 10.)

# Maximal conductances of different cell types in nS
maximal_conductances = dict(
LsoNp1a=(1000, 150, 0, 0, 1.85, 0, 7.4),
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
ninf =   (1 + exp(-(vu + 15) / 5.))**-0.5 : 1
pinf =  1. / (1 + exp(-(vu + 23) / 6.)) : 1
ntau =  ((100. / (11*exp((vu+60) / 24.) + 21*exp(-(vu+60) / 23.))) + 0.7)*ms : second
ptau = ((100. / (4*exp((vu+60) / 32.) + 5*exp(-(vu+60) / 22.))) + 5)*ms : second
"""

# Ih channel (subthreshold adaptive, non-inactivating)
eqs_ih = """
ih = ghbar*r*(Eh-v) : amp
dr/dt=q10*(rinf-r)/rtau : 1
rinf = 1. / (1+exp((vu + 76.) / 7.)) : 1
rtau = ((100000. / (237.*exp((vu+60.) / 12.) + 17.*exp(-(vu+60.) / 14.))) + 25.)*ms : second
"""

# KLT channel (low threshold K+)
eqs_klt = """
iklt = gkltbar*w**4*z*(EK-v) : amp
dw/dt=q10*(winf-w)/wtau : 1
dz/dt=q10*(zinf-z)/wtau : 1
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
ileak = gl*(El-v) : amp
"""

# h current for octopus cells
eqs_hcno = """
ihcno = gbarno*(h1*frac + h2*(1-frac))*(Eh-v) : amp
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

#eqs = """
#dv/dt = (ileak + ina + ikht + iklt + ika + ih + ihcno + I)/C : volt
#vu = v/mV : 1  # unitless v
#I = I_Bias: amp
#"""
eqs = """
dv/dt = (ileak + ina + ikht + iklt + ika + ih + ihcno + I)/C : volt
vu = v/mV : 1  # unitless v
I = I_Bias + I_Zap_Max * sin(2*pi*(m1000*(t - (t_settle+t_bias)) + f_min) * (t - (t_settle+t_bias))) : amp
#"""
#eqs = """
#dv/dt = (ileak + ina + ikht + iklt + ika + ih + ihcno + I)/C : volt
#vu = v/mV : 1  # unitless v
#if (t_settle <= t < (t_settle + sweepdur)):
#    I = I_Bias + I_Zap_Max * sin(2*pi*(m1000*(t - t_settle) + f_min) * (t - t_settle)): amp
#else:
#    I = 0: amp
#
#"""
#eqs = """
#dv/dt = (ileak + ina + ikht + iklt + ika + ih + ihcno + I)/C : volt
#vu = v/mV : 1  # unitless v
#if (t_settle <= t < (t_settle + sweepdur)):
#    I = I_Bias: amp
#else:
#    I = 0: amp
#
#"""
eqs += eqs_leak + eqs_ka + eqs_na + eqs_ih + eqs_klt + eqs_kht + eqs_hcno

#neuron = NeuronGroup(1, eqs, method='exponential_euler')
#neuron.v = El

lsonGrp = NeuronGroup(nLsons, eqs, method='exponential_euler')
# Initialize model near v_rest with no inputs
lsonGrp.v = Vrest
#vu = EL/mV # unitless v
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

#I_Gate = 1
# Current input / Zap frequency-sweep parameters
I_Bias = 0*pA
I_Zap_Max = 0*pA
sweepdur = 1*second
#sweepdur_ms = 960*ms
#f_max = 100*Hz
f_max = 2000*Hz
f_min = 1*Hz
m1000 = (f_max - f_min)/(2*(sweepdur - defaultclock.dt))
t_settle = 200*ms
t_bias = 600*ms
#zap1000_0p96 = I_Zap_Max * sin(2*pi*(m1000*(t - t_settle) + f_min) * (t - t_settle));

#M = StateMonitor(lsonGrp, 'v', record=True)
M = StateMonitor(lsonGrp, ['v','I'], record=True)
run(t_settle, report='text')  # Go to rest

I_Bias = 200.*pA
run(t_bias, report='text')  # Equilibrate to I_Bias
IBiasStr = 'IBias' + str(int(I_Bias/pA)) + 'pA'

#run(sweepdur, report='text')  # Apply Zap current
#IZapMaxStr = 'IZapMax' + str(int(I_Zap_Max/pA)) + 'pA'

I_Bias = 0.*pA
run(200*ms, report='text')

#fftStart = round((t_settle + t_bias)/defaultclock.dt)
#fftStop = round((t_settle + t_bias + sweepdur)/defaultclock.dt)
#Vfft = np.fft.fft(M[0].v[fftStart:fftStop])
#Ifft = np.fft.fft(M[0].I[fftStart:fftStop])
#Zfft = np.divide(Vfft, Ifft)

#plt.plot(M.t / ms, M[0].v / mV)
plt.plot(M.t / ms, M[0].v / mV)
plt.xlabel('t (ms)')
plt.ylabel('v (mV)')
#plt.ylim((-65,-45))
plt.show()

##plt.plot(range(5,90),abs(Vfft[5:90])) # For now, 1-Hz resolution, drop 0 Hz (large DC)
#plt.semilogx(range(10,1000),abs(Vfft[10:1000])) # For now, 1-Hz resolution, drop 0 Hz (large DC)
#plt.show()
#
##plt.plot(range(5,90),abs(Ifft[5:90])) # For now, 1-Hz resolution, drop 0 Hz (large DC)
#plt.semilogx(range(10,1000),abs(Ifft[10:1000])) # For now, 1-Hz resolution, drop 0 Hz (large DC)
#plt.show()
#
###plt.plot(range(5,90),abs(Zfft[5:90])) # For now, 1-Hz resolution, drop 0 Hz (large DC)
###plt.rcParams.update({'font.size': 18})
##plt.ylabel('impedance (MOhm)', fontsize=18)
##plt.xlabel('frequency (Hz)', fontsize=18)
##plt.semilogx(range(10,1000),abs(Zfft[10:1000])/1e6) # For now, 1-Hz resolution, drop 0 Hz (large DC)
##plt.ylim((0,200e6))
##plt.show()
#
#timeVmemFileStr = 'timeVmemNp1a_' + IBiasStr + '_' + IZapMaxStr + '.txt'
#file9 = open(timeVmemFileStr,'w')
#for index in range(len(M.t)):
#    file9.write(str(M.t[index] / ms) + " " + str(M[0].v[index] / mV) + '\n')
#file9.close()
#
#Zfft_CoreFileStr = 'ZfftNp1a_' + IBiasStr + '_' + IZapMaxStr
#absZfft_CoreFileStr = 'absZfftNp1a_' + IBiasStr + '_' + IZapMaxStr
#file10 = open(Zfft_CoreFileStr  + '.txt','w')
#file11 = open(absZfft_CoreFileStr  + '.txt','w')
#for index in range(len(Zfft)):
#    file10.write(str((index/len(Zfft)) * (1/defaultclock.dt)/Hz) + " " + str(Zfft[index]) + '\n')
#    file11.write(str((index/len(Zfft)) * (1/defaultclock.dt)/Hz) + " " + str(abs(Zfft[index])) + '\n')
#file10.close()
#file11.close()
#
#fig_absZ_Np1a, ax = plt.subplots(1, 1)
#ax.set_title('model LSO non-principal neuron (Np1a)', fontsize=16)
#ax.set_ylabel('impedance (MOhm)', fontsize=16)
#ax.set_xlabel('frequency (Hz)', fontsize=16)
#ax.set_ylim((0,200))
#ax.tick_params(axis='both', which='major', labelsize=16)
#ax.semilogx(range(10,1000),abs(Zfft[10:1000])/1e6, linewidth=2.0)
#fig_absZ_Np1a.tight_layout()
#fig_absZ_Np1a.savefig(absZfft_CoreFileStr + '.png', dpi=300)
#
##fig_out, ax = plt.subplots(1, 1)
##ax.semilogx(range(10,1000),abs(Zfft[10:1000])/1e6,  linewidth=2.0, alpha=0.3, color='gray')
##inch = 2.54
##fig_out.subplots_adjust(wspace=0.3, hspace=0.5)
##fig_out.set_size_inches(12.0/inch, 4.0 * len(electrodes)/inch)
##fig_out.tight_layout()
##fig_out.savefig(figures_path + 'ci_el_disc_fast_time_domain.pdf', dpi=300)
##fig_out.savefig(figures_path + 'ci_el_disc_fast_time_domain.png', dpi=300)
#
