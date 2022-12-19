from brian2 import *
from brian2tools import *
from brian2.units.allunits import *
import time

start_time = time.time()

#----Constants----
molar = mole/(0.001*meter**3)
mmolar = molar/1000
defaultclock.dt=.005*ms
#Area from paper
diameter = 30*umeter
length = 4.5*umeter
SchildAr = 4*pi*(diameter/2)**2
Nai = 8.9*mmolar
Nao = 145.0*mmolar
Ki = 145.0*mmolar
Ko = 5.4*mmolar
Cao = 2.0*mmolar
tauCa = 4511.0*ms
ENa = 73.0*mV
EK = -85.0*mV
gBCa = 0.000085*usiemens/SchildAr
Bi = 0.001*mmolar
ICaPmax = 0.0243*nA
INaKmax = 0.275*nA
F = 96500*coulomb/mole
T = 296*kelvin
R2 = 8.314*joule/(mole*kelvin)
DNaCa = 0.0036*mmolar**-4
KNaCa = 0.000036*nA*mmolar**-4
KMCaP = 0.0005*mmolar
KMNa = 5.46*mmolar
KMK = 0.621*mmolar
kU = 100*mmolar**-1*ms**-1
kR = 0.238*ms**-1
gamma = 0.5
r = 3
nb = 4

#----C-Type----
# gNaf = 1.95*usiemens/SchildAr
# gNas = 0.02950*usiemens/SchildAr
# gCat = 0.00035*usiemens/SchildAr
# gCan = 0.00300*usiemens/SchildAr
# gK = 0.00510*usiemens/SchildAr
# gA = 0.004*usiemens/SchildAr
# gD = 0.003*usiemens/SchildAr
# gKCa = 0.004*usiemens/SchildAr
# gBNa = 0.000525*usiemens/SchildAr
# Cm = 37.5*pF/SchildAr
# Offset Values:
# shNaf = -17.5
# shNas = -20.0
# shCat = -7.0
# shCan = -7.0
# shK = 3.0
# shA = 3.0
# shD = 3.0

#----A-Type----
gNaf = 2.05*usiemens/SchildAr
gNas = 0.00001*usiemens/SchildAr
gCat = 0.00035*usiemens/SchildAr
gCan = 0.00100*usiemens/SchildAr
gK = 0.00550*usiemens/SchildAr
gA = 0.035*usiemens/SchildAr
gD = 0.01*usiemens/SchildAr
gKCa = 0.0065*usiemens/SchildAr
gBNa = 0.000325*usiemens/SchildAr
Cm = 32.5*pF/SchildAr
#Offset Values:
shNaf = 0
shNas = 0
shCat = 0
shCan = 0
shK = 0
shA = 0
shD = 0

#----Equations----
eqs = '''
Im = -(Iion) + Istim: amp/meter**2
Iion = INaf+INas+ICat+ICan+IK+IKCa+IA+ID+(IBNa+IBCa)+ICaP+INaK-2.*INaCa: amp/meter**2
Istim: amp/meter**2
dv/dt = (Im)/Cm : volt

#Unitless Vm
vu = v/mV: 1
'''

eqs_Cat = '''
ICat = gCat*dtt*ft*(v-ECa): amp/meter**2

###gating variable dt renamed to dtt to avoid conflict with internal variable dt###
ddtt/dt = (dttinf - dtt)/taudtt: 1
dft/dt = (ftinf - ft)/tauft: 1

taudtt = (22.0*exp(-(0.052)**2*(vu+68.0+shCat)**2)+2.5)*ms: second
tauft = (103.0*exp(-(0.050)**2*(vu+58.0+shCat)**2)+12.5)*ms: second

dttinf = 1.0/(1.0+exp((vu+54.00+shCat)/-5.75)): 1
ftinf = 1.0/(1.0+exp((vu+68.00+shCat)/6.0)): 1
'''

eqs_Can = '''
ICan = gCan*dn*(0.55*fn1 + 0.45*fn2)*(v-ECa): amp/meter**2

ddn/dt = (dninf - dn)/taudn: 1
dfn1/dt = (fn1inf - fn1)/taufn1: 1
dfn2/dt = (fn2inf - fn2)/taufn2: 1

taudn = (3.25*exp(-(0.042**2)*(vu+31.0+shCan)**2)+0.395)*ms: second
taufn1 = (33.5*exp(-(0.0395**2)*(vu+30.0+shCan)**2) + 5.0)*ms: second
taufn2 = (225.0*exp(-(0.0275**2)*(vu+40.0+shCan)**2)+75.0)*ms: second

dninf = 1.0/(1.0+exp((vu+20.0+shCan)/-4.5)): 1
fn1inf = 1.0/(1.0+exp((vu+20.0+shCan)/25.0)): 1
fn2inf = rn + 1.0/(1.0+exp((vu+40+shCan)/10.0)): 1
rn = 0.2/(1.0+exp((vu+5.0+shCan)/-10.0)): 1
'''

eqs_Naf = '''
INaf = gNaf*mf**3*hf*j*(v-ENa): amp/meter**2

dmf/dt = (mfinf-mf)/taumf: 1
dhf/dt = (hfinf-hf)/tauhf: 1
dj/dt = (jinf-j)/tauj: 1

taumf = (0.75*exp(-(0.0635**2)*(vu+40.35+shNaf)**2)+0.12)*ms: second
tauhf = (6.5*exp(-(0.0295**2)*(vu+75.00+shNaf)**2)+0.55)*ms: second
tauj = (25.0/(1.0+exp((vu-20.00+shNaf)/4.5))+0.01)*ms: second

mfinf = 1.0/(1.0+exp((vu+41.35+shNaf)/-4.75)): 1
hfinf = 1.0/(1.0+exp((vu+62.00+shNaf)/4.50)): 1
jinf = 1.0/(1.0+exp((vu+40.00)/1.50)): 1
'''

eqs_Nas = '''
INas = gNas*mss**3*hs*(v-ENa): amp/meter**2

dmss/dt = (mssinf-mss)/taumss: 1
dhs/dt = (hsinf - hs)/tauhs: 1


taumss = (1.50*exp(-(0.0595**2)*(vu+20.35+shNas)**2)+0.15)*ms: second
tauhs = (4.95*exp(-(0.0335)**2*(vu+20.00+shNas)**2)+0.75)*ms: second


mssinf = 1.0/(1.0+exp((vu+20.35+shNas)/-4.45)): 1
hsinf = 1.0/(1.0+exp((vu+18.00+shNas)/4.50)): 1


'''

eqs_K = '''
IK = gK*n*(v-EK): amp/meter**2
dn/dt = (ninf-n)/taun: 1

taun = (1.0/(an+bn))+1.0*ms: second

ninf = 1.0/(1.0+exp((vu+14.62+shK)/-18.38)): 1

an = (0.001265*(vu+14.273)/(1.0-exp((vu+14.273+shK)/-10.0)))/ms: Hz
bn = 0.125*exp((vu+55.0+shK)/-2.5)/ms: Hz
'''

eqs_A = '''
IA = gA*p**3*q*(v-EK): amp/meter**2

dp/dt = (pinf-p)/taup: 1
dq/dt = (qinf-q)/tauq: 1

taup = (5.0*exp(-(0.022**2)*(vu+65.0+shA)**2)+2.5)*ms: second
tauq = (100.0*exp(-(0.035**2)*(vu+30.0+shA)**2)+10.5)*ms: second

pinf = 1.0/(1.0+exp((vu+28.0+shA)/-28.0)): 1
qinf = 1.0/(1.0+exp((vu+58.0+shA)/7.0)): 1
'''

eqs_D = '''
###gating variables x and y renamed to xx and yy to avoid conflict with internal variables x and y
ID = gD*xx**3*yy*(v-EK): amp/meter**2

dxx/dt = (xxinf-xx)/tauxx: 1
dyy/dt = (yyinf-yy)/tauyy: 1

tauxx = (5.0*exp(-(0.022)**2*(vu+65.0)**2)+2.5)*ms: second
tauyy = 7500.0*ms: second

xxinf = 1.0/(1.0+exp((vu+39.59+shD)/-14.68)): 1
yyinf = 1.0/(1.0+exp((vu+48.0+shD)/7.0)): 1
'''

eqs_KCa = '''
IKCa = gKCa*c*(v-EK): amp/meter**2

dc/dt = (cinf-c)/tauc: 1

tauc = 4.5/(ac+bc): second

cinf = ac/(ac+bc): 1

ac = (750.0*(Cai/mmolar)*exp((vu-10.)/12.0))/ms: Hz
bc = (0.05*exp((vu-10)/-60.0))/ms: Hz
'''

eqs_BNa = '''
IBNa = gBNa*(v-ENa): amp/meter**2
'''

eqs_BCa = '''
IBCa = gBCa*(v-ECa): amp/meter**2
'''

eqs_NaCa = '''
INaCa = KNaCa*(DFin-DFout)/S/surfAr: amp/meter**2
surfAr = 4*pi*(diameter/2)**2: meter**2
S = 1.0+DNaCa*(Cai*Nao**r+Cas*Nai**r): 1
DFin = Nai**r*Cas*exp((r-2)*gamma*v*F/(R2*T)): mole**4/meter**12
DFout = Nao**r*Cai*exp((r-2)*(gamma-1)*v*F/(R2*T)): mole**4/meter**12
'''

eqs_CaP = '''
ICaP = ICaPmax * (Cai/(Cai+KMCaP))/surfAr: amp/meter**2
'''

eqs_NaK = '''
INaK = INaKmax * (Nai/(Nai+KMNa))**3*(Ko/(Ko+KMK))**2*((vu+150)/(vu+200))/surfAr: amp/meter**2
'''

eqs_CaConc = '''
#internal volume
#10% of cell volume for organelles
voli = (4/3*pi*(diameter/2)**3)*0.9: meter**3

#perineural space volume, thickness 1.0 um (pg 2345) 
#vols = length*pi*((diameter/2)+1*um)**2-voli: meter**3
vols = 4/3*pi*(diameter/2+1*um)**3-(voli/0.9): meter**3

dOC/dt = kU*Cai*(1-OC)-kR*OC: 1
diffOC = kU*Cai*(1-OC)-kR*OC: Hz

#diffOC = dOC/dt in Hz for dCai/dt equation
dCai/dt = (2*INaCa-ICan-ICat-IBCa-ICaP)*surfAr/(2*voli*F) - nb*Bi*diffOC: mole/meter**3
dCas/dt = (Cas - Cao)/tauCa + (-2*INaCa+ICan+ICat+IBCa+ICaP)*surfAr/(2*vols*F): mole/meter**3

ECa = (R2*T/(2*F))*log(Cas/Cai)-78.7*mV: volt
'''


eqs = eqs+eqs_Cat+eqs_Can+eqs_Naf+eqs_Nas+eqs_K+eqs_A+eqs_D+eqs_KCa+eqs_BNa+eqs_BCa+eqs_NaCa+eqs_CaP+eqs_NaK+eqs_CaConc


#----Simulation----
neuron = NeuronGroup(1, eqs,clock=Clock(defaultclock.dt, name='clock_1'),
                    threshold='v > -40*mV',
                    refractory='v > -40*mV', name = 'neuron')
                    
#Restore state from holding voltage steady state and release holding voltage
restore(filename='/Users/Edward/Documents/Ind Study/schild-73')

# monitor1=StateMonitor(neuron, ('v','IA','ID','IKCa','IK','ICan', 'INaf', 'Cas', 'Cai', 'INaK', 'INaCa', 'IBNa', 'IBCa'), record=True, name = 'monitor1')
# monitor2=StateMonitor(neuron, ('dtt', 'ft', 'dn', 'fn1', 'fn2', 'mf', 'hf', 'j', 'mss', 'hs', 'n', 'p', 'q', 'xx', 'yy', 'c', 'Cai', 'Cas', 'OC'), record=True, name = 'monitor2')

#Stimulus
neuron.Istim[:] = (40*pA)/SchildAr
run(200*ms, report='text')

#Turn off stimulus
neuron.Istim[:] = 0*uA/cm**2
run(100*ms, report='text')

print("Finishd in %s seconds" % (time.time() - start_time))


show()