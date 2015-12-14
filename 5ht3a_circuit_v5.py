

import csv
import neuron
from neuron import h
from neuron import gui
import matplotlib.pyplot as plt
import subprocess
from neuron import load_mechanisms
import random
import math
import numpy as np
import itertools
# just adding some comments
class cell(object):
	def __init__(self):

		self.make_sections()
		self.connsect()
		self.biophys()

	def make_sections(self):
		self.soma=neuron.h.Section(name='soma',cell=self)
		self.dend=neuron.h.Section(name='dend',cell=self)
	def connsect(self):
		self.dend.connect(self.soma(0.5))
		
	def biophys(self):
		
		self.soma.nseg=1
		self.soma.diam=10
		self.soma.L=10
		self.soma.Ra=125
		self.soma.insert("hh")
		

      		#self.soma.cascale_aabBK = .01

		
		#self.soma.ca_i_Ca=.1
		self.dend.nseg=5
		self.dend.L=300
		self.dend.diam=0.5
		self.dend.Ra=125
		self.dend.insert("pas")
		self.soma.gkbar_hh=.038

		







def set_stim(seg,amp,delay,dur):
	stim=neuron.h.IClamp(seg(0.5))
	stim.delay=delay
	stim.dur=dur
	stim.amp=amp
	return stim

def recvolt(seg):
	rec_v=neuron.h.Vector()
	rec_v.record(seg._ref_v)
	return rec_v

def rectime(seg):
	rec_t=neuron.h.Vector()
	rec_t.record(seg._ref_t)
	return rec_t

load_mechanisms('/home/ben/Desktop/5ht3a_project/neuron_model/downloaded_neuron_models/Hipp_paper_code/')

h.xopen('/home/ben/Desktop/5ht3a_project/neuron_model/downloaded_neuron_models/Hipp_paper_code/basket_cell17S.hoc')
h.xopen('/home/ben/Desktop/5ht3a_project/neuron_model/downloaded_neuron_models/Hipp_paper_code/pyramidal_cell_14Vb.hoc')
h.xopen('/home/ben/Desktop/5ht3a_project/neuron_model/downloaded_neuron_models/Hipp_paper_code/olm_cell2.hoc')


pyrstren=np.arange(.005,.01,.001)
pvstren=np.arange(0,0.3,.05)


olmexcstren=np.arange(.05,.055,.001)
olmstren=np.arange(.04,.05,.002)
olmcombs=itertools.product(olmexcstren,olmstren)
pvpyrcombs=itertools.product(pyrstren,pvstren)

numiter=100
ppp=1
for ggg in olmcombs:
	print 'run#'+str(ppp)+'/30'
	ppp=ppp+1


	spikecount=np.zeros((numiter,3))

	for hhh in range(0,numiter):

		sccell=cell()
		rec_v_sccell=recvolt(sccell.soma(0.5))
		sccell_stim=set_stim(sccell.soma,1,100,2)
		rec_v_sccell=recvolt(sccell.soma(0.5))

		pyrcell=h.PyramidalCell()
		pyrcell.soma.Ra=30
		#pyrstim=set_stim(pyrcell.soma,2)

		rec_v_pyrcell=recvolt(pyrcell.soma(0.5))

		synpyr=neuron.h.Exp2Syn(pyrcell.radTprox(0))
		synpyr.tau1=0.5
		synpyr.tau2=3
		synpyr.e=0

		pvin=h.BasketCell()
		pvin.soma.Ra=200
		pvin.soma.cm=1
		#pvstim=set_stim(pvin.soma,2)


		rec_v_pvin=recvolt(pvin.soma(0.5))
		synpv=neuron.h.Exp2Syn(pvin.radT2(0))
		synpvpyr=neuron.h.Exp2Syn(pyrcell.soma(0.5))
		synpv.tau1=0.5
		synpv.tau2=3
		synpv.e=0
		synpvpyr.tau1=0.5
		synpvpyr.tau2=3
		synpvpyr.e=-75
		###SPECIFY THE VALUES OF THE SECTION TO BE RECORDED##
		#record time from NEURON (neuron.h._ref_t)
		rec_t=h.Vector()
		rec_t.record(neuron.h._ref_t)


		spontpyrepsc=neuron.h.Exp2Syn(pyrcell.soma(0.5))
		spontpyrepsc.tau1=0.5
		spontpyrepsc.tau2=3
		spontpyrepsc.e=0
		spontpyripsc=neuron.h.Exp2Syn(pyrcell.soma(0.5))
		spontpyripsc.tau1=1
		spontpyripsc.tau2=8
		spontpyripsc.e=-75



		#######Inititate OLM cell and define synapses##############
		olmcell=h.OLMCell()
		rec_v_olm=recvolt(olmcell.soma(0.5))
		#olmstim=set_stim(olmcell.soma,.1)

		synolm=neuron.h.Exp2Syn(olmcell.dend1(0))
		synolm.tau1=0.5
		synolm.tau2=3
		synolm.e=0

		synolmpv=neuron.h.Exp2Syn(pvin.soma(0.5))
		synolmpv.tau1=2
		synolmpv.tau2=10
		synolmpv.e=-75

		synolmpyr=neuron.h.Exp2Syn(pyrcell.soma(0.5))
		synolmpyr.tau1=2
		synolmpyr.tau2=10
		synolmpyr.e=-75

		spontolmepsc=neuron.h.Exp2Syn(olmcell.soma(0.5))
		spontolmepsc.tau1=0.5
		spontolmepsc.tau2=3
		spontolmepsc.e=0
		spontolmipsc=neuron.h.Exp2Syn(olmcell.soma(0.5))
		spontolmipsc.tau1=1
		spontolmipsc.tau2=8
		spontolmipsc.e=-75
		#############################################################



		#####################synaptic connections########################################
		####lag between EPSP and IPSP should be ~3 ms from Karnup and Stelzer (1999, J. Physiol)
		ampasyn=neuron.h.NetCon(sccell.soma(0.5)._ref_v,synpyr,0,0,.008,sec=sccell.soma)#original conductance=.006; somewhere between .0042 and .0045 pyr cells fire 100%
		#probability of exciting pyr with no inhibition: .0044=.94, .0045=.98, .0046=.97,.0047=.97,.005=.95,.006=.98,.0065=.99,.007=1,.01=1
		ampasynpv=neuron.h.NetCon(sccell.soma(0.5)._ref_v,synpv,0,0,0.04,sec=sccell.soma)#original weight=.03
		gabasynpvpyr=neuron.h.NetCon(pvin.soma(0.5)._ref_v,synpvpyr,-40,1,.2,sec=pvin.soma) #original weight=.01
		ampasynolm=neuron.h.NetCon(sccell.soma(0.5)._ref_v,synolm,-20,0,ggg[0],sec=sccell.soma)#original weight=.01
		gabasynolmpv=neuron.h.NetCon(olmcell.soma(0.5)._ref_v,synolmpv,-30,0,ggg[1],sec=olmcell.soma)#original weight=.04
		gabasynolmpyr=neuron.h.NetCon(olmcell.soma(0.5)._ref_v,synolmpyr,-30,0,0,sec=olmcell.soma)#original weight=.01
		###############################################################

		###spontaneous events on pyramidal cell######
		#generate poisson process giving sepsc and sipsc times in time course
		rateconst = float(2000)
		timekeep=0
		events=[]
		eventtimespyrepsc=[]
		for i in range(1,14003):
			ppval=1.0-math.exp(-(i-timekeep)/rateconst)
			if ppval>random.random():
				events.append(1)
				eventtimespyrepsc.append(i)
				timekeep=i
			else:
				events.append(0)
				timekeep=i
		spontpyrepsc=list()
		for jj in range(0,len(eventtimespyrepsc)):
			spontpyrepsc.append(neuron.h.AlphaSynapse(pyrcell.soma(0.5)))
			spontpyrepsc[jj].onset=(eventtimespyrepsc[jj]/14002.0)*350.0
			#print "eventimes[jj]="+str(spontpyrepsc[jj].onset)
			spontpyrepsc[jj].e=0
			spontpyrepsc[jj].gmax=.002*np.random.lognormal(.002,.2,1)
			spontpyrepsc[jj].tau=1


		timekeep=0
		events=[]
		eventtimespyripsc=[]
		for i in range(1,14003):
			ppval=1.0-math.exp(-(i-timekeep)/rateconst)
			if ppval>random.random():
				events.append(1)
				eventtimespyripsc.append(i)
				timekeep=i
			else:
				events.append(0)
				timekeep=i
		spontpyripsc=list()
		for jj in range(0,len(eventtimespyripsc)):
			spontpyripsc.append(neuron.h.AlphaSynapse(pyrcell.soma(0.5)))
			spontpyripsc[jj].onset=(eventtimespyripsc[jj]/14002.0)*350.0
			#print "eventimes[jj]="+str(spontpyripsc[jj].onset)
			spontpyripsc[jj].e=-75
			spontpyripsc[jj].gmax=.002*np.random.lognormal(.002,.2,1)
			spontpyripsc[jj].tau=5

		######################################################

		####Spontaneous synaptic events on PV cell#############

		timekeep=0
		events=[]
		eventtimespvepsc=[]
		for i in range(1,14003):
			ppval=1.0-math.exp(-(i-timekeep)/rateconst)
			if ppval>random.random():
				events.append(1)
				eventtimespvepsc.append(i)
				timekeep=i
			else:
				events.append(0)
				timekeep=i
		spontpvepsc=list()
		for jj in range(0,len(eventtimespvepsc)):
			spontpvepsc.append(neuron.h.AlphaSynapse(pvin.soma(0.5)))
			spontpvepsc[jj].onset=(eventtimespvepsc[jj]/14002.0)*350.0
			#print "eventimes[jj]="+str(spontpvepsc[jj].onset)
			spontpvepsc[jj].e=0
			spontpvepsc[jj].gmax=.002*np.random.lognormal(.002,.2,1)
			spontpvepsc[jj].tau=1

		rateconst=float(2000)
		timekeep=0
		events=[]
		eventtimespvipsc=[]
		for i in range(1,14003):
			ppval=1.0-math.exp(-(i-timekeep)/rateconst)
			if ppval>random.random():
				events.append(1)
				eventtimespvipsc.append(i)
				timekeep=i
			else:
				events.append(0)
				timekeep=i
		spontpvipsc=list()
		for jj in range(0,len(eventtimespvipsc)):
			spontpvipsc.append(neuron.h.AlphaSynapse(pvin.soma(0.5)))
			spontpvipsc[jj].onset=(eventtimespvipsc[jj]/14002.0)*350.0
			#print "eventimes[jj]="+str(spontpvipsc[jj].onset)
			spontpvipsc[jj].e=-75
			spontpvipsc[jj].gmax=.002*np.random.lognormal(.002,.2,1)
			spontpvipsc[jj].tau=5


		######################################################

		####Spontaneous synaptic events on OLM cell#############

		timekeep=0
		events=[]
		eventtimesolmepsc=[]
		for i in range(1,14003):
			ppval=1.0-math.exp(-(i-timekeep)/rateconst)
			if ppval>random.random():
				events.append(1)
				eventtimesolmepsc.append(i)
				timekeep=i
			else:
				events.append(0)
				timekeep=i
		spontolmepsc=list()
		for jj in range(0,len(eventtimesolmepsc)):
			spontolmepsc.append(neuron.h.AlphaSynapse(olmcell.soma(0.5)))
			spontolmepsc[jj].onset=(eventtimesolmepsc[jj]/14002.0)*350.0
			#print "eventimes[jj]="+str(spontolmepsc[jj].onset)
			spontolmepsc[jj].e=0
			spontolmepsc[jj].gmax=.002*np.random.lognormal(.002,.2,1)
			spontolmepsc[jj].tau=1

		rateconst=float(2000)
		timekeep=0
		events=[]
		eventtimesolmipsc=[]
		for i in range(1,14003):
			ppval=1.0-math.exp(-(i-timekeep)/rateconst)
			if ppval>random.random():
				events.append(1)
				eventtimesolmipsc.append(i)
				timekeep=i
			else:
				events.append(0)
				timekeep=i
		spontolmipsc=list()
		for jj in range(0,len(eventtimesolmipsc)):
			spontolmipsc.append(neuron.h.AlphaSynapse(olmcell.soma(0.5)))
			spontolmipsc[jj].onset=(eventtimesolmipsc[jj]/14002.0)*350.0
			#print "eventimes[jj]="+str(spontolmipsc[jj].onset)
			spontolmipsc[jj].e=-75
			spontolmipsc[jj].gmax=.002*np.random.lognormal(.002,.2,1)
			spontolmipsc[jj].tau=5



		######################################################

		######constant currents############
		#set_stim(seg,amp,delay,duration)
		pyrstim=set_stim(pyrcell.soma,0,0,0)
		pvstm=set_stim(pvin.soma,0,0,0)
		olmcellstim=set_stim(olmcell.soma,-.5,0,350)

		#################################

		h.finitialize(-65)
		#set time of the simulation
		tstop=350

		#Run the simulation
		neuron.run(tstop)

		time=rec_t.as_numpy()

		pvin_volt=rec_v_pvin.as_numpy()
		pyrcell_volt=rec_v_pyrcell.as_numpy()
		olmcell_volt=rec_v_olm.as_numpy()
		sccell_volt=rec_v_sccell.as_numpy()
		'''
		plt.plot(time,olmcell_volt,color='b',linewidth=3.0)
		plt.show()
	
		print "length pvin_volt="+str(len(pvin_volt))
	
		plt.plot(time,pvin_volt,color='r',linewidth=3.0)
		#plt.plot(time,olmcell_volt,color='b',linewidth=3.0)
		plt.plot(time,sccell_volt,color='g',linewidth=3.0)
		plt.plot(time,pyrcell_volt,color='k',linewidth=3.0)
		plt.show()
		'''


		'''
		print "pv Rm="+str(pvin.soma.Ra)
		print "pyr Rm="+str(pyrcell.soma.Ra)
		print "length of v vectors:"+str(len(pyrcell_volt))
		print "mean of events="+str(np.mean(events))
		print "length of events="+str(len(events))


		print "length of spont epsc="+str(len(spontpyrepsc))
		'''
		'''
		print "max value pyr:"+str(max(pyrcell_volt[4081:4400]))#search for max in (4081:4400)*350/14002=102~110 ms
		print "max value pv:"+str(max(pvin_volt[4081:4400]))
		print "max value pv:"+str(max(olmcell_volt[4081:4400]))
		print "iteration:"+str(hhh+1)
		'''
		if max(pyrcell_volt[4081:4400])>-20:
			spikecount[hhh,0]=1
		if max(pvin_volt[4081:4400])>0:
			spikecount[hhh,1]=1
		if max(olmcell_volt[4081:4400])>0:
			spikecount[hhh,2]=1
	'''
	print "pyrstren="+str(ggg[0])
	print "pvstren="+str(ggg[1])
	'''
	print "mean pyr:"+str(np.mean(spikecount[:,0]))
	
	with open('test2.csv','a') as testfile:
		csvwriter=csv.writer(testfile)
		csvwriter.writerow([ggg[0],ggg[1],np.mean(spikecount[:,0])])

	
	print "mean pv:"+str(np.mean(spikecount[:,1]))
	print "mean olm:"+str(np.mean(spikecount[:,2]))



