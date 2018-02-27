# XRD pattern calculation for Heusler alloys and related materials
# Copyright (C) 2018  Patrick R. LeClair
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.


#for L21 X2YZ compound -or- C1b XYZ -or- DO19 X2YZ
#or (X1X2)(Y1Y2)(Z1Z2) to have quaternaries, etc
#L21 = space group no. 225, Fm-3m
#C1b = space group no. 216, F-43m
# ^ they have the same general rules allowed hkl, only sites change a bit
#DO19 = space group no. 194, P63/mmc
#very different from other two ... 

import csv
from numpy import *
import numpy as np
import matplotlib.pyplot as plt

#class so we can do things like Site.a2, Site.h6, etc
class positions: 
	pass
	
######## HERE IS A BLOCK THAT YOU EDIT ########

A = 5.784 #5.22 							#a lattice constant angstroms
C = 4.24
hmax = kmax = lmax = 4
THETA_MAX = 120
THETA_MIN = 5 
XRAY = "Co"   # or "Cu"

if (XRAY=="Co"): #Co				#specify what x-ray wavelength to use
	Lambda = 1.79					#in Angstroms
	print("Co Ka")
elif (XRAY=="Cu"): #Cu
	Lambda = 1.54
	print("Cu Ka")
else: 
	Lambda = 1.54	
	print("No X-ray wavelength specified, defaulting to Cu Ka.")

#which structure?
structure = "DO19" #L21, C1b, or DO19 are coded

if (structure=="L21"):
	C1b=0
	DO19=0
	print("L21 structure")
elif (structure=="C1b"):
	L21=0
	DO19=0
	print("C1b structure")
elif (structure=="DO19"):
	L21=0
	C1b=0
	print("DO19 structure")
	

x = 1.0/6 #x,z parameters for structure, if needed
z = 1.0/6

######## STOP EDITING HERE AND SCROLL DOWN ########

#up to three elements present. can share a site by assigning occupancy 
#e.g., X = [Fe,Sites.h6,2.0/3] and Y = [Mn,Sites.h6,1.0/3] to have a 2:1 Fe-Mn mix on 6h

X1 = []
X2 = []
Y1 = []
Y2 = []
Z1 = []
Z2 = []

#Wyckoff positions
Sites = positions()

#naming convention: wycoff letter + multiplicity b/c we can't use names like Site.2a
if (structure=="L21"):
	Sites.a4 = ['a4', (0.0,0.0,0.0)]
	Sites.b4 = ['b4', (0.5,0.5,0.5)]
	Sites.c8 = ['c8', (1.0/4,1.0/4,1.0/4), (1.0/4,1.0/4,3.0/4)]
	Sites.d24 = ['d24', (0.0,0.25,0.25), (0.0,0.75,0.25), (0.25,0,0,0.25), (0.25,0.0,0.75), (0.25,0.25,0.0), (0.75,0.25,0.0)]
	SitesTuple = [Sites.a4,Sites.b4,Sites.c8,Sites.d24]

if (structure=="C1b"):
	Sites.a4 = ['a4', (0.0,0.0,0.0)]
	Sites.b4 = ['b4', (0.5,0.5,0.5)]
	Sites.c4 = ['c4', (0.25,0.25,0.25)]
	Sites.d4 = ['d4', (0.75,0.75,0.75)]
	SitesTuple = [Sites.a4,Sites.b4,Sites.c4,Sites.d4]
	
if (structure=="DO19"):
	Sites.a2 = ['a2', (0.0,0.0,0.0), (0.0,0.0,0.5)]
	Sites.b2 = ['b2', (0.0,0.0,0.25), (0.0,0.0,0.75)]
	Sites.c2 = ['c2', (1.0/3,2.0/3,0.25), (2.0/3,1.0/3,0.75)]
	Sites.d2 = ['d2', (1.0/3,2.0/3,0.75), (2.0/3,1.0/3,0.25)]
	Sites.e4 = ['e4', (0.0,0.0,z), (0.0,0.0,z+0.5), (0.0,0.0,-z), (0,0,0,-z+0.5)]
	Sites.f4 = ['f4', (1.0/3,2.0/3,z), (2.0/3,1.0/3,z+0.5), (2.0/3,1.0/3,-z), (1.0/3,2.0/3,-z+0.5)]
	Sites.g6 = ['g6', (0.5,0.0,0.0), (0.0,0.5,0.0), (0.5,0.5,0.0), (0.5,0.0,0.5), (0.0,0.5,0.5), (0.5,0.5,0.5)]
	Sites.h6 = ['h6', (x,2*x,0.25), (-2*x,-x,0.25), (x,-x,0.25), (-x,-2*x,0.75), (2*x,x,0.75), (-x,x,0.75)]
	SitesTuple = [Sites.a2,Sites.b2,Sites.c2,Sites.d2,Sites.e4,Sites.f4,Sites.g6,Sites.h6]

#element data: atomic scattering factors, atomic numbers
Al = 13  #these assignments are to make the scattering factor matrices readable
Si = 14
Ti = 22
V  = 23
Cr = 24
Mn = 25
Fe = 26
Co = 27
Ni = 28
Cu = 29
Zn = 30
Ga = 31
Ge = 32
Sn = 50
Sb = 51
Gd = 64

elements = {      		#this structure is to make output readable mostly
	'13':'Al',
	'14':'Si',
	'22':'Ti',
	'23':'V',
	'24':'Cr',
	'25':'Mn',
	'26':'Fe',
	'27':'Co',
	'28':'Ni',
	'29':'Cu',
	'30':'Zn',
	'31':'Ga',
	'32':'Ge',
	'50':'Sn',
	'51':'Sb',
	'64':'Gd'
	}

#5th order polynomial to approximate fo vs sin(T)/labmda=1/2d for elements 0-92
#array element [i][6] is f', element [i][7] is f" (dispersion corrections)
#use http://it.iucr.org/Cb/ch6o1v0001/sec6o1o1/ table 6.1.1.1 for fo 
#use http://it.iucr.org/Cb/ch4o2v0001/sec4o2o6/ table 4.2.6.8 for f', f"

#TODO CHECK numbers below again against ITC once more

#do NOT average structure factors to mix elements. 
#since intensity depends on f^2, this will not work right. instead, adjust occupancy.
#can have e.g. X = [Fe,Sites.h6,2.0/3] and Y = [Mn,Sites.h6,1.0/3] to randomize sites

ScatteringFactor = [[0 for x in range(8)] for x in range(111)]  #roentgenium

ScatteringFactor[13][0] = 13.000			#Al
ScatteringFactor[13][1] = 1.5744
ScatteringFactor[13][2] = -347.93
ScatteringFactor[13][3] = 1974.7
ScatteringFactor[13][4] = -4533.9
ScatteringFactor[13][5] = 3760.6
if (XRAY=="Co"):
	ScatteringFactor[13][6] = 0.2551
	ScatteringFactor[13][7] = 0.3276*1j
elif (XRAY=="Cu"):
	ScatteringFactor[13][6] = 0.2130
	ScatteringFactor[13][7] = 0.2455*1j
else:#default to Cu Ka if XRAY ill-defined
	ScatteringFactor[13][6] = 0.2130
	ScatteringFactor[13][7] = 0.2455*1j	
	
ScatteringFactor[14][0] = 13.986			#Si
ScatteringFactor[14][1] = 3.2804
ScatteringFactor[14][2] = -375.64
ScatteringFactor[14][3] = 1954.7
ScatteringFactor[14][4] = -4126.7
ScatteringFactor[14][5] = 3179.2
if (XRAY=="Co"):
	ScatteringFactor[14][6] = 0.2979
	ScatteringFactor[14][7] = 0.4384*1j
elif (XRAY=="Cu"):
	ScatteringFactor[14][6] = 0.2541
	ScatteringFactor[14][7] = 0.3302*1j
else:#default to Cu Ka if XRAY ill-defined
	ScatteringFactor[14][6] = 0.2541
	ScatteringFactor[14][7] = 0.3302*1j		

ScatteringFactor[22][0] = 22.027			#Ti
ScatteringFactor[22][1] = -1.4931
ScatteringFactor[22][2] = -423.94
ScatteringFactor[22][3] = 2286
ScatteringFactor[22][4] = -5328.8
ScatteringFactor[22][5] = 4628.8
if (XRAY=="Co"):
	ScatteringFactor[22][6] = -0.0617
	ScatteringFactor[22][7] = 2.3213*1j
elif (XRAY=="Cu"):
	ScatteringFactor[22][6] = 0.2191
	ScatteringFactor[22][7] = 1.8069*1j
else:
	ScatteringFactor[22][6] = 0.2191
	ScatteringFactor[22][7] = 1.8069*1j	
	
ScatteringFactor[23][0] = 23.018			#V
ScatteringFactor[23][1] = -0.3458
ScatteringFactor[23][2] = -422.51
ScatteringFactor[23][3] = 2194.2
ScatteringFactor[23][4] = -4992.1
ScatteringFactor[23][5] = 4262.4
if (XRAY=="Co"):
	ScatteringFactor[23][6] = -0.3871
	ScatteringFactor[23][7] = 2.6994*1j
elif (XRAY=="Cu"):
	ScatteringFactor[23][6] = 0.0687
	ScatteringFactor[23][7] = 2.1097*1j
else:
	ScatteringFactor[23][6] = 0.0687
	ScatteringFactor[23][7] = 2.1097*1j	
	
ScatteringFactor[24][0] = 24.012			#Cr
ScatteringFactor[24][1] = -0.2731
ScatteringFactor[24][2] = -340.71
ScatteringFactor[24][3] = 1512.3
ScatteringFactor[24][4] = -3120
ScatteringFactor[24][5] = 2533.2
if (XRAY=="Co"):
	ScatteringFactor[24][6] = -0.9524
	ScatteringFactor[24][7] = 3.1130*1j
elif (XRAY=="Cu"):
	ScatteringFactor[24][6] = -0.1635
	ScatteringFactor[24][7] = 2.4439*1j
else:
	ScatteringFactor[24][6] = -0.1635
	ScatteringFactor[24][7] = 2.4439*1j		

ScatteringFactor[25][0] = 25.005			#Mn
ScatteringFactor[25][1] = 1.0076
ScatteringFactor[25][2] = -405.86
ScatteringFactor[25][3] = 1968.7
ScatteringFactor[25][4] = -4290.2
ScatteringFactor[25][5] = 3555.6
if (XRAY=="Co"):
	ScatteringFactor[25][6] = -2.0793
	ScatteringFactor[25][7] = 3.5546*1j
elif (XRAY=="Cu"):	
	ScatteringFactor[25][6] = -0.5299
	ScatteringFactor[25][7] = 2.8052*1j
else:	
	ScatteringFactor[25][6] = -0.5299
	ScatteringFactor[25][7] = 2.8052*1j	

ScatteringFactor[26][0] = 26.001 			#Fe
ScatteringFactor[26][1] = 1.4161
ScatteringFactor[26][2] = -394.85
ScatteringFactor[26][3] = 1856.6
ScatteringFactor[26][4] = -3968.7
ScatteringFactor[26][5] = 3246.3
if (XRAY=="Co"):
	ScatteringFactor[26][6] = -3.3307
	ScatteringFactor[26][7] = 0.4901*1j
elif (XRAY=="Cu"):
	ScatteringFactor[26][6] = -1.1336
	ScatteringFactor[26][7] = 3.1974*1j
else: #default to Cu Ka if undefined
	ScatteringFactor[26][6] = -1.1336
	ScatteringFactor[26][7] = 3.1974*1j

ScatteringFactor[27][0] = 26.998 			#Co
ScatteringFactor[27][1] = 1.6884
ScatteringFactor[27][2] = -382.6
ScatteringFactor[27][3] = 1744.9
ScatteringFactor[27][4] = -3658.7
ScatteringFactor[27][5] = 2953.2
if (XRAY=="Co"):
	ScatteringFactor[27][6] = -2.023
	ScatteringFactor[27][7] = 0.5731*1j
elif (XRAY=="Cu"):
	ScatteringFactor[27][6] = -2.3653
	ScatteringFactor[27][7] = 3.6143*1j
else: #default to Cu Ka if undefined
	ScatteringFactor[27][6] = -2.3653
	ScatteringFactor[27][7] = 3.6143*1j	

ScatteringFactor[28][0] = 27.996			#Ni
ScatteringFactor[28][1] = 1.8709
ScatteringFactor[28][2] = -370.04
ScatteringFactor[28][3] = 1638.4
ScatteringFactor[28][4] = -3371.7
ScatteringFactor[28][5] = 2865.9
if (XRAY=="Co"):
	ScatteringFactor[28][6] = -1.5664
	ScatteringFactor[28][7] = 0.6662*1j
elif (XRAY=="Cu"):
	ScatteringFactor[28][6] = -3.0029
	ScatteringFactor[28][7] = 0.5091*1j
else: #default to Cu Ka if undefined
	ScatteringFactor[28][6] = -3.0029
	ScatteringFactor[28][7] = 0.5091*1j		
	
ScatteringFactor[29][0] = 29.000			#Cu
ScatteringFactor[29][1] = 0.8695
ScatteringFactor[29][2] = -290.18
ScatteringFactor[29][3] = 1092.2
ScatteringFactor[29][4] = -2046.2
ScatteringFactor[29][5] = 1570.7
if (XRAY=="Co"):
	ScatteringFactor[29][6] = -1.2789
	ScatteringFactor[29][7] = 0.7700*1j
elif (XRAY=="Cu"):
	ScatteringFactor[29][6] = -1.9646
	ScatteringFactor[29][7] = 0.5888*1j
else: #default to Cu Ka if undefined
	ScatteringFactor[29][6] = -1.9646
	ScatteringFactor[29][7] = 0.5888*1j			

ScatteringFactor[30][0] = 29.993			#Zn
ScatteringFactor[30][1] = 2.0638
ScatteringFactor[30][2] = -345.01
ScatteringFactor[30][3] = 1442.4
ScatteringFactor[30][4] = -2860.7
ScatteringFactor[30][5] = 2220
if (XRAY=="Co"):
	ScatteringFactor[30][6] = -1.0843
	ScatteringFactor[30][7] = 0.8857*1j
elif (XRAY=="Cu"):	
	ScatteringFactor[30][6] = -1.5491
	ScatteringFactor[30][7] = 0.6778*1j
else:	
	ScatteringFactor[30][6] = -1.5491
	ScatteringFactor[30][7] = 0.6778*1j		

ScatteringFactor[31][0] = 30.997			#Ga
ScatteringFactor[31][1] = 1.6692
ScatteringFactor[31][2] = -391.71
ScatteringFactor[31][3] = 1769.8
ScatteringFactor[31][4] = -3634.4
ScatteringFactor[31][5] = 2838.8
if (XRAY=="Co"):
	ScatteringFactor[31][6] = -0.9200
	ScatteringFactor[31][7] = 1.0138*1j
elif (XRAY=="Cu"):	
	ScatteringFactor[31][6] = -1.2846
	ScatteringFactor[31][7] = 0.7736*1j
else:	
	ScatteringFactor[31][6] = -1.2846
	ScatteringFactor[31][7] = 0.7736*1j	

ScatteringFactor[32][0] = 31.987			#Ge
ScatteringFactor[32][1] = 3.2051
ScatteringFactor[32][2] = -441.19
ScatteringFactor[32][3] = 2014.7
ScatteringFactor[32][4] = -4086
ScatteringFactor[32][5] = 3128.6
if (XRAY=="Co"):
	ScatteringFactor[32][6] = -0.7781
	ScatteringFactor[32][7] = 1.557*1j
elif (XRAY=="Cu"):	
	ScatteringFactor[32][6] = -1.0885
	ScatteringFactor[32][7] = 0.8855*1j
else:	
	ScatteringFactor[32][6] = -1.0885
	ScatteringFactor[32][7] = 0.8855*1j	
	
ScatteringFactor[50][0] = 49.989			#Sn
ScatteringFactor[50][1] = 3.5716
ScatteringFactor[50][2] = -619.41
ScatteringFactor[50][3] = 2713.6
ScatteringFactor[50][4] = -5432.4
ScatteringFactor[50][5] = 4210.2
if (XRAY=="Co"):
	ScatteringFactor[50][6] = -0.3097
	ScatteringFactor[50][7] = 6.9896*1j
elif (XRAY=="Cu"):
	ScatteringFactor[50][6] = 0.0259
	ScatteringFactor[50][7] = 5.4591*1j
else:
	ScatteringFactor[50][6] = 0.0259
	ScatteringFactor[50][7] =  5.4591*1j

ScatteringFactor[51][0] = 50.975			#Sb
ScatteringFactor[51][1] = 5.443
ScatteringFactor[51][2] = -664.91
ScatteringFactor[51][3] = 2895.8
ScatteringFactor[51][4] = -5644.3
ScatteringFactor[51][5] = 4213.4
if (XRAY=="Co"):
	ScatteringFactor[51][6] = -0.5189
	ScatteringFactor[51][7] = 7.5367*1j
elif (XRAY=="Cu"):
	ScatteringFactor[51][6] = -0.0562
	ScatteringFactor[51][7] = 5.8946*1j
else:
	ScatteringFactor[51][6] = -0.0562
	ScatteringFactor[51][7] = 5.8946*1j	
	
ScatteringFactor[64][0] = 64.058			#Gd
ScatteringFactor[64][1] = -5.0071
ScatteringFactor[64][2] = -653.45
ScatteringFactor[64][3] = 3073.1
ScatteringFactor[64][4] = -6751.1
ScatteringFactor[64][5] = 5752.2
if (XRAY=="Co"):
	ScatteringFactor[64][6] = -9.3863
	ScatteringFactor[64][7] = 3.9016*1j
elif (XRAY=="Cu"):
	ScatteringFactor[64][6] = -8.8380
	ScatteringFactor[64][7] = 11.9157*1j
else:
	ScatteringFactor[64][6] = -8.8380
	ScatteringFactor[64][7] = 11.9157*1j

#generate atomic scattering factor with the data in the matrices above
#we are not doing the Debye-Waller factor. or the thickness/absorption correction
#Debye-Waller data is hard to come by, and we presume bulk samples here.

def f(element,d):  			
	s = 1.0/(2.0*d)    #sin(theta)/lambda
	f=0
	for i in range(6):
		f += (ScatteringFactor[element][i])*(s**i)  #poly approx of non-dispersive part
	f += ScatteringFactor[element][6] + ScatteringFactor[element][7] #add dispersion
	return (f)		#note element [7] is imaginary, so return value is complex

#general rules for a given space group on allowed hkl
			
def rules(h,k,l):	#general rules for allowed hkl 
	allowed=1
	if (structure=="L21" or structure=="C1b"): #(225 and 216 have same rules)
		if (((h+k)%2==1) or ((h+l)%2==1) or ((l+k)%2==1)):
			allowed = 0
		if ((h==0) and ((k%2==1) or (l%2==1))):
			allowed=0
		if ((h==k) and ((h+l)%2==1)):
			allowed=0
		if ((k==0 and l==0) and (h%2==1)):
			allowed=0
		if (h==0 and k==0 and l==0):
			allowed=0
	if (structure=="DO19"):
		if ((h==k) and (l%2 == 1)):
			allowed=0
		if ((h==0 and k==0) and (l%2 == 1)):
			allowed=0
		if (h==0 and k==0 and l==0):
			allowed=0
	return(allowed)

#calculate structure factor, including rules for specific sites in a space group

def F_hkl(site,h,k,l):		
	F=0
	if (structure=="L21" and (site[0]=='c8' or site[0]=='d24')): 
		if (h%2==0):
			for i in range (1,len(site)):
				F += exp(2*pi*1j*(site[i][0]*h+site[i][1]*k+site[i][2]*l))
	elif (structure=="DO19" and (site[0]=='a2' or site[0]=='b2' or site[0]=='e4' or site[0]=='g6')):
		if (l%2==0):
			for i in range (1,len(site)):
				F += exp(2*pi*1j*(site[i][0]*h+site[i][1]*k+site[i][2]*l))
	elif (structure=="DO19" and (site[0]=='c2' or site[0]=='d2' or site[0]=='f4')):
		if ((((l%2==0) or ((h-k)%3==1) or ((h-k)%3==2)))):
			for i in range (1,len(site)):
				F += exp(2*pi*1j*(site[i][0]*h+site[i][1]*k+site[i][2]*l))			
	else: #any other site
		for i in range (1,len(site)):
				F += exp(2*pi*1j*(site[i][0]*h+site[i][1]*k+site[i][2]*l))
	return (F)

#find d spacing, bragg angle, and Lorentz-polarization factors

def d_hkl(h,k,l):								#cubic d spacing for hkl
	if (structure=="L21" or structure=="C1b"):
		tmp = (h**2+k**2+l**2)/(A*A)
		return (1.0/sqrt(tmp))
	elif (structure=="DO19"):
		tmp = (4.0/3.0)*(h*h+h*k+k*k)/(A*A) + l*l/(C*C)
		return (1.0/sqrt(tmp))
	
def bragg(d,Lambda):
	tmp = (Lambda/(2.0*d))
	if tmp<=1:
		angle = 2.0*degrees(arcsin(tmp)) 
	else: 
		angle = 0	#bad arcsin
	return(angle)
	
def Lorentz_Pol(d,Lambda):
	if ((Lambda/(2.0*d)) <= 1.0):
		theta = arcsin(Lambda/(2.0*d)) #radians!
		LP = (1.0+(cos(2.0*theta))**2)/(sin(theta)*sin(2.0*theta))	
	else:
		theta = 0	#bad arcsin
		LP = 0
	return(LP)

#find all the peaks. no need for multiplicity factor since we brute force all combos
			
def Pattern(X1,X2,Y1,Y2,Z1,Z2,plot,outputfile,outputsites):			
	h=-hmax
	k=-kmax
	l=-lmax

	F_X1=0
	F_Y1=0
	F_Z1=0
	F_X2=0
	F_Y2=0
	F_Z2=0

	pattern = []
	allowed = 0

	while (h<=hmax):
		while (k<=kmax):
			while (l<=lmax):
				i = -(h+k) 	#fourth hexagonal index
				allowed=rules(h,k,l)
				if (allowed):
					#e.g. 50% occupied, scale accordingly
					#X has X[element,site,occupancy]
					F_X1 = F_hkl(X1[1],h,k,l)*X1[2] 
					F_X2 = F_hkl(X2[1],h,k,l)*X2[2]	
					F_Y1 = F_hkl(Y1[1],h,k,l)*Y1[2] 
					F_Y2 = F_hkl(Y2[1],h,k,l)*Y2[2]	
					F_Z1 = F_hkl(Z1[1],h,k,l)*Z1[2] 
					F_Z2 = F_hkl(Z2[1],h,k,l)*Z2[2]
					d = d_hkl(h,k,l)
					two_theta = bragg(d,Lambda)
					LP = Lorentz_Pol(d,Lambda)
					I = LP*(absolute(F_X1*f(X1[0],d)+F_Y1*f(Y1[0],d)+F_Z1*f(Z1[0],d)+F_X2*f(X2[0],d)+F_Y2*f(Y2[0],d)+F_Z2*f(Z2[0],d)))**2
					if ((two_theta<THETA_MAX) and (two_theta>THETA_MIN)):
						pattern.append([two_theta,h,k,i,l,F_X1,F_X2,F_Y1,F_Y2,F_Z1,F_Z2,I,d])
				allowed=0
				l+=1
			k+=1
			l=-lmax
		h+=1
		k=-kmax

	pattern.sort(key=lambda x: x[0])	#sort list on 2-theta value

	if (outputfile):
		OutFile = "./output/"+elements[str(X1[0])]+X1[1][0]+elements[str(X2[0])]+X2[1][0]+elements[str(Y1[0])]+Y1[1][0]+elements[str(Y2[0])]+Y2[1][0]+elements[str(Z1[0])]+Z1[1][0]+elements[str(Z2[0])]+Z2[1][0]+"scattering-factors-hkil"+".csv"
		of = open(OutFile, 'wt')
		DataOut=csv.writer(of)
		
		DataOut.writerow(["Elements"])
		DataOut.writerow(["X1 X2 Y1 Y2 Z1 Z2"])
		DataOut.writerow([(elements[str(X1[0])],elements[str(X2[0])],elements[str(Y1[0])],elements[str(Y2[0])],elements[str(Z1[0])],elements[str(Z2[0])])])		
		DataOut.writerow(["Sites"])
		DataOut.writerow(["X1 X2 Y1 Y2 Z1 Z2"])
		DataOut.writerow([X1[1][0],X2[1][0],Y1[1][0],Y2[1][0],Z1[1][0],Z2[1][0]])
		DataOut.writerow(["Occupancy"])
		DataOut.writerow(["X1 X2 Y1 Y2 Z1 Z2"]) 
		DataOut.writerow([X1[2],X2[2],Y1[2],Y2[2],Z1[2],Z2[2]])
	
		DataOut.writerow(['2T','h','k','i','l','Fx1','Fy1','Fz1','Fx2','Fy2','Fz2','I','d'])
		DataOut.writerows(pattern)
		of.close()    	

	if (outputlistverbose):
		print "data for all allowed (hkl)"
		if (structure=="DO19"):
			print "\n\n2Theta \t hkil \t X1 \t X2 \t Y1 \t Y2 \t Z1 \t Z2 \t I \t d (A)"
			for x in pattern:
				print '{0:8.2f}\t ({1:1d},{2:1d},{3:1d},{4:1d}) \t {5:6.2f} \t {6:6.2f}  \t {7:8.2f} \t {8:6.2f} \t {9:6.2f}  \t {10:8.2f} \t {11}  \t {12:8.3f} '.format(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9],x[10],x[11],x[12])
		else:
			print "\n\n2Theta \t hkl \t X1 \t X2 \t Y1 \t Y2 \t Z1 \t Z2 \t I \t d (A)"
			for x in pattern:
				print '{0:8.2f}\t ({1:1d},{2:1d},{3:1d}) \t {4:6.2f} \t {5:6.2f}  \t {6:8.2f} \t {7:6.2f} \t {8:6.2f}  \t {9:8.2f} \t {10}  \t {11:8.3f} '.format(x[0], x[1], x[2], x[3], x[5], x[6], x[7], x[8], x[9],x[10],x[11],x[12])

	pattern_dict = {}
	for t in pattern:
		if t[0] in pattern_dict:
			pattern_dict[t[0]] = pattern_dict[t[0]]+t[11]
		else:
			pattern_dict[t[0]] = t[11]

# 	could uncomment this to just print intensity in AU and skip normalization
#	print "\n2Theta \t I (approx)" #raw intensity, not normalized
#	for key,value in sorted(pattern_dict.items()):
#		print '{0:8.2f} {1:8.2f}'.format(key,value)		
#		DataOut.writerow([key, value])
 
#	print "\n"

	# This section scales the intensities relative to the maximum and prints the results

	pattern_dict2 = pattern_dict  #creates new dictionary so I don't write over old one
 
	lenint = len(sorted(pattern_dict2.items())) #determines number of peaks

	maxi = 0  #sets the initial maximum peak value to 0
	maxtheta = 0

	for i in range(lenint):           			#determines the largest peak 
		if sorted(pattern_dict2.items())[i][1]>=maxi:
			maxi = sorted(pattern_dict2.items())[i][1]
			maxtheta = sorted(pattern_dict2.items())[i][0]

	for i in range(lenint):    					#scales all of the peaks 
		pattern_dict2[sorted(pattern_dict2.items())[i][0]] = 100*pattern_dict2[sorted(pattern_dict2.items())[i][0]]/maxi

	if (outputlist):
		print "\n2Theta \t I (normalized)"      #prints the results

		for key,value in sorted(pattern_dict2.items()):
			print '{0:8.2f}\t{1:8.6f}'.format(key,value)		

	if (outputfile):
		OutFile = "./output/"+elements[str(X1[0])]+X1[1][0]+elements[str(X2[0])]+X2[1][0]+elements[str(Y1[0])]+Y1[1][0]+elements[str(Y2[0])]+Y2[1][0]+elements[str(Z1[0])]+Z1[1][0]+elements[str(Z2[0])]+Z2[1][0]+"peak-list"+".csv"
		of = open(OutFile, 'at')
		DataOut=csv.writer(of)
		DataOut.writerow(['2T','I (norm.)'])
		for key,value in sorted(pattern_dict2.items()):
			DataOut.writerow([key, value])
		of.close()    

	if (outputsites):	
		print "\nElements X1=%s X2=%s Y1=%s Y2=%s Z1=%s Z2=%s"%(elements[str(X1[0])],elements[str(X2[0])],elements[str(Y1[0])],elements[str(Y2[0])],elements[str(Z1[0])],elements[str(Z2[0])])		
		print "Sites X1=%s X2=%s Y1=%s Y2=%s Z1=%s Z2=%s \nOccupancy X1=%s X2=%s Y1=%s Y2=%s Z1=%s Z2=%s \n "%(X1[1][0],X2[1][0],Y1[1][0],Y2[1][0],Z1[1][0],Z2[1][0],X1[2],X2[2],Y1[2],Y2[2],Z1[2],Z2[2])	
		print "Max peak at %f\n"%maxtheta	

 	if (plot or plotfile):						#prepare data to plot
		x = [] 
		y = []
		for key,value in sorted(pattern_dict2.items()):
			x.append(key)
			y.append(value)
	
	if (plot):									#popup plot
		plt.bar(x,y,width=0.2,color='b',edgecolor='b')
		plt.axis([15, 120, 0, 100])
		plt.show()
	if (plotfile):								#just dump a PNG
		plt.bar(x,y,width=0.2,color='b',edgecolor='b')
		plt.axis([15, 120, 0, 100])
		foo = "./output/"+elements[str(X1[0])]+X1[1][0]+elements[str(X2[0])]+X2[1][0]+elements[str(Y1[0])]+Y1[1][0]+elements[str(Y2[0])]+Y2[1][0]+elements[str(Z1[0])]+Z1[1][0]+elements[str(Z2[0])]+Z2[1][0]+"graph"+".png"
		plt.savefig(foo)
	return(maxtheta)
	
	
######## HERE IS A BLOCK THAT YOU EDIT ########	

#various output options
plot=1						#pop-up plot of xrd pattern
plotfile=0					#save a png of the Pattern
outputfile=1				#output scattering factors & peak list to file
outputlist=1				#print a summary list of peaks to the tty
outputlistverbose=1			#print the ENTIRE list of peaks to the tty
outputsites=1				#print which elements are on which sites to tty
search_xyz=0				#search over x/y/z parameter 
search_sites=0				#try switching site assignments [can take a while]

#X has [element,site,occupancy]
#if e.g., X1 and X2 elements share a site, take care that total occupancy isn't > 1

# X1 = [Fe,Sites.b4,1.0]    		#1 atom			#set as C1b structure above!
# X2 = [Co,Sites.d4,1.0] 			#1 atom
# Y1 = [Co,Sites.c4,0.5]			#1 atom
# Y2 = [Ti,Sites.c4,0.5]			#1 atom
# Z1 = [Ge,Sites.a4,0.5]			#1 atom
# Z2 = [Ge,Sites.a4,0.5]			#1 atom

X1 = [Fe,Sites.h6,2.0/3]   	#set as DO19 structure above! define x & z above!
X2 = [Fe,Sites.h6,0] 
Y1 = [Mn,Sites.h6,1.0/3]
Y2 = [Mn,Sites.h6,0]
Z1 = [Ge,Sites.d2,1.0]
Z2 = [Ge,Sites.d2,0]

Pattern(X1,X2,Y1,Y2,Z1,Z2,plot,outputfile,outputsites)		





	
######## BELOW THIS YOU ARE IN THE WEEDS ########
######## needs to be modified to handle multiple atoms on each posn, e.g., X1 X2 etc#####

# CODE BELOW to search brute for varying x parameter on 6h site in DO19

# if (search_xyz):
# 	for j in range (0,7):
# 		x = j/6.0
# 		print "z=%f"%(x)
# 		Sites.h6 = ['h6', (x,2*x,1.0/4), (-2*x,-x,1.0/4), (x,-x,1.0/4), (-x,-2*x,3.0/4), (2*x,x,3.0/4), (-x,x,3.0/4)]
# 		X = [Fe,Sites.h6,2.0/3]
# 		Y = [Mn,Sites.h6,1.0/3]
# 		Z = [Ge,Sites.c2,1.0]
# 		maxtheta=Pattern(X,Y,Z,plot,outputfile,outputsites)	
# 		if (maxtheta>30.0):
# 			print "maxtheta %f for X=%s Y=%s Z=%s and x=%f"%(maxtheta,X[1][0],Y[1][0],Z[1][0],x)


### Code snippets below do various brute force searches. you are in the weeds now.

#Pattern(X,Y,Z,plot,outputfile,outputsites)

#search_sites=1

# if (search_sites):
# 	for i in range (len(SitesTuple)):
# 		X = [FeCo,SitesTuple[i],num_X/(len(SitesTuple[i])-1)]
# 		for j in range (len(SitesTuple)):
# 			if (j!=i):
# 				Y = [Ge,SitesTuple[j],num_Y/(len(SitesTuple[j])-1)]
# 			for k in range (len(SitesTuple)):
# 				if (j!=i and k!=j and k!=i):
# 					Z = [Ti,SitesTuple[k],num_Z/(len(SitesTuple[k])-1)]
# 					maxtheta=Pattern(X,Y,Z,plot,outputfile,outputsites)	
# 					if (maxtheta>30.0):
# 						print "maxtheta %f for X=%s Y=%s Z=%s"%(maxtheta,X[1][0],Y[1][0],Z[1][0])
# 						#print "occ %f %f %f"%(num_X/(len(SitesTuple[i])-1),num_Y/(len(SitesTuple[j])-1),num_Z/(len(SitesTuple[k])-1))

# CODE BELOW to search brute force over xyz for specific peaks
# search_xyz=0
# 
# if (search_xyz):
# 	for k in range (0,6):
# 		y = k/10.0
# 		for i in range (0,6): 	
# 			z = i/10.0	
# 			for j in range (0,6):
# 				x = j/10.0
# 				Sites.f4 = ['f4', (1.0/3,2.0/3,z), (2.0/3,1.0/3,z+0.5), (2.0/3,1.0/3,-z), (1.0/3,2.0/3,-z+0.5)]
# 				Sites.g6 = ['g6', (0.5,0.0,0.0), (0.0,0.5,0.0), (0.5, 0.5, 0.0), (0.5,0.0,0.5), (0.0,0.5,0.5), (0.5,0.5,0.5)]
# 				Sites.h6 = ['h6', (x,2*x,1.0/4), (-2*x,-x,1.0/4), (x,-x,1.0/4), (-x,-2*x,3.0/4), (2*x,x,3.0/4), (-x,x,3.0/4)]
# 				Sites.j12 = ['j12', (x,y,1.0/4), (-y,x-y,1.0/4), (-x+y,-x,1.0/4), (-x,-y,3.0/4), (y,-x+y,3.0/4), (x-y,x,3.0/4), (y,x,3.0/4), (x-y,-y,3.0/4), (-x,-x+y,3.0/4), (-y,-x,1.0/4), (-x+y,y,1.0/4), (x,x-y,1.0/4)]
# 				X = [FeCo,Sites.h6,1.0/3]
# 				Y = [Ti,Sites.f4,1.0/2]
# 				Z = [Ge,Sites.g6,1.0/3]
# 				maxtheta=Pattern(X,Y,Z,plot,outputfile,outputsites)	
# 				if (maxtheta>50.0):
# 					print "maxtheta %f for X=%s Y=%s Z=%s"%(maxtheta,X[1][0],Y[1][0],Z[1][0])		
# 					print "x=%f, y=%f, z=%f"%(x,y,z)

#num_X = 2.0
#num_Y = 2.0
#num_Z = 2.0
	
#Pattern(X,Y,Z,plot,outputfile,outputsites)

#search_sites=1

# if (search_sites):
# 	for i in range (len(SitesTuple)):
# 		X = [FeCo,SitesTuple[i],num_X/(len(SitesTuple[i])-1)]
# 		for j in range (len(SitesTuple)):
# 			if (j!=i):
# 				Y = [Ge,SitesTuple[j],num_Y/(len(SitesTuple[j])-1)]
# 			for k in range (len(SitesTuple)):
# 				if (j!=i and k!=j and k!=i):
# 					Z = [Ti,SitesTuple[k],num_Z/(len(SitesTuple[k])-1)]
# 					maxtheta=Pattern(X,Y,Z,plot,outputfile,outputsites)	
# 					if (maxtheta>30.0):
# 						print "maxtheta %f for X=%s Y=%s Z=%s"%(maxtheta,X[1][0],Y[1][0],Z[1][0])
# 						#print "occ %f %f %f"%(num_X/(len(SitesTuple[i])-1),num_Y/(len(SitesTuple[j])-1),num_Z/(len(SitesTuple[k])-1))
