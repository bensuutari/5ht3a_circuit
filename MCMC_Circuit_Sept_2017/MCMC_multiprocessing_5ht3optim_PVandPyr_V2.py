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
from sim_circuit_multiprocessing import *
startnewsim=True#start simulation from scratch=True, continue simulating from previous data=False
numiter=20
ppp=1
MCsimnum=1000
params=np.zeros((MCsimnum,2))
numiters=20
if startnewsim:
	params[0,0]=0.016843152324863357 #0.002 0.014751983753824395,0.081784536364045771
	params[0,1]=0.074018763154778931 #0.016 
else:
	with open('5ht3optim_PVandPyr_BuildHist_v2.csv','rb') as testfile:
		csvreader=csv.reader(testfile)
		for row in csvreader:
			prevdata.append(row)
	params[0,0]=float(prevdata[-1][1]) #orignal .04
	params[0,1]=float(prevdata[-1][2])	#original .2
	prevdata=list()
	print prevdata[-1][1]
	print prevdata[-1][2]

####Get the template cells from Cutsuridis et al., 2010#######
load_mechanisms('Hipp_paper_code/')

h.xopen('Hipp_paper_code/basket_cell17S.hoc')
h.xopen('Hipp_paper_code/pyramidal_cell_14Vb.hoc')
h.xopen('Hipp_paper_code/olm_cell2.hoc')
##########################################################




#variables to store spike probability
meanpyr=np.zeros((MCsimnum,1))
meanpv=np.zeros((MCsimnum,1))

#run first simulation to make initial guess
spikecount=np.zeros((numiters,2))

paramslist=[(params[0,0],params[0,1]) for yyy in range(0,numiters)]
pool=multiprocessing.Pool()
dummyvar=pool.map(runsims,paramslist)

dummypyr=[dummyvar[i][0] for i in range(0,numiters)]
dummypv=[dummyvar[i][1] for i in range(0,numiters)]
pool.close()
meanpyr[0]=np.mean(dummypyr)
meanpv[0]=np.mean(dummypv)
print "mean pyr:"+str(meanpyr[0])



print "mean pv:"+str(meanpv[0])
#print "mean olm:"+str(np.mean(spikecount[:,2]))
#order for datadump: SCpyr, SCPV, PVpyr, Pyrprob,PVprob

with open('5ht3optim.csv','a') as testfile:
	csvwriter=csv.writer(testfile)
	csvwriter.writerow([.0065,float(params[1,0]),float(params[1,0]),float(meanpyr[0]),float(meanpv[0])])





#use MCMC to generate pdf

for ggg in range(1,MCsimnum):
	print 'entered loop'
	print 'MCMC sim#'+str(ggg)
	#print 'run#'+str(ppp)+'/'+str(len(olmexcstren)*len(olmstren))
	ppp=ppp+1


	spikecount=np.zeros((numiters,2))
	#step the parameters
	#first get some random numbers
	step5ht3pyr=random.random()/1000
	step5ht3pv=random.random()/100
	if random.random()>0.5:
		step5ht3pyr=-step5ht3pyr
	if random.random()>0.5:
		step5ht3pv=-step5ht3pv
	
	if (params[ggg-1,0]+step5ht3pyr)<=0:
		step5ht3pyr=-step5ht3pyr
	if (params[ggg-1,1]+step5ht3pv)<=0:
		step5ht3pv=-step5ht3pv
	
	if (params[ggg-1,0]+step5ht3pyr)>.002*10:
		step5ht3pyr=-step5ht3pyr
	if (params[ggg-1,1]+step5ht3pv)>.016*10:
		step5ht3pv=-step5ht3pv
	paramsc=[params[ggg-1,0]+step5ht3pyr,params[ggg-1,1]+step5ht3pv]
	 
	paramslist=[paramsc for yyy in range(0,numiters)]
	pool=multiprocessing.Pool()
	dummyvar=pool.map(runsims,paramslist)

	dummypyr=[dummyvar[i][0] for i in range(0,numiters)]
	dummypv=[dummyvar[i][1] for i in range(0,numiters)]
	dummyvar=None
	meanpyr[ggg]=np.mean(dummypyr)
	meanpv[ggg]=np.mean(dummypv)
	pool.close()


	print "mean pyr:"+str(meanpyr[ggg])
	print "mean pv:"+str(meanpv[ggg])
	#print "mean olm:"+str(np.mean(spikecount[:,2]))
	print('Cost='+str((np.exp(-4*abs((meanpyr[ggg]-0.8)))+np.exp(-4*abs((meanpv[ggg]-0.4))))/2))
	if random.random()<(np.exp(-4*abs((meanpyr[ggg]-0.8)))+np.exp(-4*abs((meanpv[ggg]-0.4))))/2:
		params[ggg,0]=paramsc[0]
		params[ggg,1]=paramsc[1]
	else:
		params[ggg,0]=params[ggg-1,0]
		params[ggg,1]=params[ggg-1,1]
		meanpyr[ggg]=meanpyr[ggg-1]
		meanpv[ggg]=meanpv[ggg-1]
		

	print params[0:ggg,:]
	with open('5ht3optim_PVandPyr_BuildHist_v2.csv','a') as testfile:
		csvwriter=csv.writer(testfile)
		csvwriter.writerow([.0065,params[ggg,0],params[ggg,1],float(meanpyr[ggg]),float(meanpv[ggg])])



