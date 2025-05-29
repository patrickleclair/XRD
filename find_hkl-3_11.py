# XRD pattern calculation for Heusler alloys and related materials
# Copyright (C) 2018-2023  Patrick R. LeClair
# this version should work with python 3
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
#or (X1X2)(Y1Y2)(Z1Z2) to have quaternaries, etc, or Xa inverse Heusler
#obviously also handles rocksalt, fcc
#C1b = space group no. 216, Fm-3m. structures B3 C1b C15b
#L21 = space group no. 225, F-43m. structures fcc, B1, C1, DO3, D2f, D84, L21, Ca7Ge
#see https://homepage.univie.ac.at/michael.leitner/lattice/spcgrp/cubic.html
# ^ they have the same general rules allowed hkl, only sites change a bit
#DO19 = space group no. 194, P63/mmc
#covers hcp, graphite, and many other structures
#https://homepage.univie.ac.at/michael.leitner/lattice/spcgrp/hexagonal.html
#DO22 = space group no. 139, I4/mmm. structures DO22, DO23, A6 (monotomic L1_0)
#see https://homepage.univie.ac.at/michael.leitner/lattice/spcgrp/tetragonal.html
#space group 46 = orthorhombic, adopted by FeTiSi 

#todo: automatic calculation of linear absorption coefficient mu for thickness dep.
	#need another mu parameter for each element & x-ray used
	#and code to do the composition weighting
	#for now: global MU we assume is proper avg you want at the x-ray of interest
#todo: parameters for other elements
#todo: other space groups, variation of parameters associated
#todo: DW factor is a hack, make a separate function for this for clarity?

import csv
from numpy import *
import numpy as np
import matplotlib.pyplot as plt

#class so we can do things like Site.a2, Site.h6, etc
class positions: 
	pass
	
######## HERE IS A BLOCK THAT YOU EDIT ########

#specimen types to adjust lorentz-polarization factor		
POWDER = 0 #don't change this, read on
SINGLE_XTAL = 1	

#calculation ranges
hmax = kmax = lmax = 10
THETA_MAX = 120
THETA_MIN = 5 

#wavelength and corrections
XRAY = "Cu"   	# "Co" or "Cu"; no others implemented currently
DISPERSION = 1  #include f' and f" dispersion corrections to atomic scattering factor?
				#we turn this off to compare with other software
				
DEBYE_WALLER = 1 		#include debye-waller correction or no 
SAMPLE_TYPE = SINGLE_XTAL 	# POWDER or SINGLE_XTAL, to determine Lorentz-polarization 
FILM = 1
THICKNESS = 19.5e-7 #in cm

MU = 3031  #in 1/cm, for thickness corr.

#lattice constants
A = 6.00 #5.784 #5.22 					#a lattice constant angstroms
B = 10.97
C = 6.37 #4.24

#which space group to generate structure?
space_group = "SG225" #SG46 SG224, SG216, SG194, SG139

###############################################

if (SAMPLE_TYPE == POWDER):
	print("Powder Lorentz-polarization correction")
elif (SAMPLE_TYPE == SINGLE_XTAL):
	print("Single crystal Lorentz-polarization correction")
	
if (DISPERSION):
	print("Dispersion corrections to atomic scattering factor f' and f'' included")
else:
	print("Dispersion corrections to atomic scattering factor f' and f'' NOT included")

if (DEBYE_WALLER):
	print("Debye-Waller corrections included")
else:
	print("Debye-Waller corrections NOT included")
	
if (FILM):
	print("Thin film: thickness correction applied, t=%s cm, mu=%s 1/cm"%(THICKNESS,MU))
else:
	print("Bulk assumed, no thickness correction")

if (XRAY=="Co"): #Co				#specify what x-ray wavelength to use
	Lambda = 1.79026				#in Angstroms
	print("Co Ka")					#TODO: may as well include other wavelengths
elif (XRAY=="Cu"): #Cu
	Lambda = 1.54184
	print("Cu Ka")
else: 
	Lambda = 1.54184	
	print("No X-ray wavelength specified, defaulting to Cu Ka.")	

print("%s structure"%(space_group))
print("a lattice parameter %s A"%(A))
if (space_group=="SG46"):
	print("b lattice parameter %s A"%(B))
if (space_group=="SG194" or space_group=="SG139"):
	print("c lattice parameter %s A"%(C))

x = 0.0 #x,z parameters for structure, if needed
y = 1.0/6
z = 0.12

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
if (space_group=="SG46"):
	Sites.c8 = ['c8', (x,y,z), (-x,y,z), (x+0.5,-y,z), (-x+0.5,y,z)]
	Sites.b4 = ['b4', (0.25,y,z), (0.75,-y,z)]
	Sites.a4 = ['a4', (0,0,z), (0.5,0,z)]

if (space_group=="SG225"):
	Sites.a4 = ['a4', (0.0,0.0,0.0)]
	Sites.b4 = ['b4', (0.5,0.5,0.5)]
	Sites.c8 = ['c8', (1.0/4,1.0/4,1.0/4), (1.0/4,1.0/4,3.0/4)]
	Sites.d24 = ['d24', (0.0,0.25,0.25), (0.0,0.75,0.25), (0.25,0,0,0.25), (0.25,0.0,0.75), (0.25,0.25,0.0), (0.75,0.25,0.0)]
	SitesTuple = [Sites.a4,Sites.b4,Sites.c8,Sites.d24]

if (space_group=="SG224"): #2nd origin choice
	Sites.a2 = ['a2', (0.0,0.0,0.0), (0.5,0.5,0.5)]
	Sites.b4 = ['b4', (1.0/4,1.0/4,1.0/4), (3.0/4,3.0/4,1.0/4), (3.0/4,1.0/4,3.0/4), (1.0/4,3.0/4,3.0/4)]
	Sites.c4 = ['c4', (3.0/4,3.0/4,3.0/4), (1.0/4,1.0/4,3.0/4), (1.0/4,3.0/4,1.0/4), (3.0/4,1.0/4,1.0/4)]
	Sites.d6 = ['d6', (0,0.5,0.5), (0.5,0,0.5), (0.5,0.5,0), (0,0.5,0), (0.5,0,0), (0,0,0.5)]
	SitesTuple = [Sites.a2,Sites.b4,Sites.c4,Sites.d6]

if (space_group=="SG216"):
	Sites.a4 = ['a4', (0.0,0.0,0.0)]
	Sites.b4 = ['b4', (0.5,0.5,0.5)]
	Sites.c4 = ['c4', (0.25,0.25,0.25)]
	Sites.d4 = ['d4', (0.75,0.75,0.75)]
	SitesTuple = [Sites.a4,Sites.b4,Sites.c4,Sites.d4]
	
if (space_group=="SG194"):
	Sites.a2 = ['a2', (0.0,0.0,0.0), (0.0,0.0,0.5)]
	Sites.b2 = ['b2', (0.0,0.0,0.25), (0.0,0.0,0.75)]
	Sites.c2 = ['c2', (1.0/3,2.0/3,0.25), (2.0/3,1.0/3,0.75)]
	Sites.d2 = ['d2', (1.0/3,2.0/3,0.75), (2.0/3,1.0/3,0.25)]
	Sites.e4 = ['e4', (0.0,0.0,z), (0.0,0.0,z+0.5), (0.0,0.0,-z), (0,0,0,-z+0.5)]
	Sites.f4 = ['f4', (1.0/3,2.0/3,z), (2.0/3,1.0/3,z+0.5), (2.0/3,1.0/3,-z), (1.0/3,2.0/3,-z+0.5)]
	Sites.g6 = ['g6', (0.5,0.0,0.0), (0.0,0.5,0.0), (0.5,0.5,0.0), (0.5,0.0,0.5), (0.0,0.5,0.5), (0.5,0.5,0.5)]
	Sites.h6 = ['h6', (x,2*x,0.25), (-2*x,-x,0.25), (x,-x,0.25), (-x,-2*x,0.75), (2*x,x,0.75), (-x,x,0.75)]
	SitesTuple = [Sites.a2,Sites.b2,Sites.c2,Sites.d2,Sites.e4,Sites.f4,Sites.g6,Sites.h6]
	
if (space_group=="SG139"):
	Sites.a2 = ['a2', (0.0,0.0,0.0)]
	Sites.b2 = ['b2', (0.0,0.0,0.5)]
	Sites.c4 = ['c4', (0.0,0.5,0.0), (0.5,0.0,0.0)]
	Sites.d4 = ['d4', (0.0,0.5,0.25), (0.5,0.0,0.25)]
	Sites.e4 = ['e4', (0.0,0.0,z), (0.0,0.0,-z)]
	SitesTuple = [Sites.a2,Sites.b2,Sites.c4,Sites.d4,Sites.e4]
	
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
FeCo = 100

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
	'64':'Gd',
	'100':'FeCo',
	}

# 5 gaussian approximation to f(s) from ref below. this is what VESTA does fyi.
# Acta Cryst. (1995). A51,416-431
# New Analytical Scattering-Factor Functions for Free Atoms and Ions
# BY D. WAASMAIER AND A. KIRFEL
# a[i] are elements 0-4, b[i] are elements 5-9, c is element 10
# f = sum a[i] * exp[-b[i]*s^2] + c where s=sin(theta)/lambda 
# alternatively see data at 
# http://it.iucr.org/Cb/ch6o1v0001/sec6o1o1/ table 6.1.1.1 for fo 
# and approximate as you will
# Note Fe and Co are tough b/c near absorption for Co, Cu Ka radiation. Careful!

# dispersion corrections are:
# array element [i][11] is f', element [i][12] is f" (dispersion corrections)
# linear interpolation of https://physics.nist.gov/PhysRefData/FFast/html/form.html
# INCLUDES nuclear thompson and relativistic corrections to f1
# so (my f') = f' + f_NT = f1 + f_rel + f_NT - Z as they put it
# can also use http://it.iucr.org/Cb/ch4o2v0001/sec4o2o6/ table 4.2.6.8 for f', f"
# ends up being very close

# do NOT average structure factors to mix elements. 
# since intensity depends on f^2, this will not work right. instead, adjust occupancy.
# can have e.g. X = [Fe,Sites.h6,2.0/3] and Y = [Mn,Sites.h6,1.0/3] to randomize sites

# last item in list [13] is Debye-Waller factor B in (angstrom)^2 from 
# DebyeWaller Factors and Absorptive Scattering Factors of Elemental Crystals
# L.-M. Peng, G. Ren, S. L. Dudarev and M. J. Whelan
# Acta Crystallographica Section A Foundations of Crystallography
# or otherwise where noted
# we use data from elemental crystals (not ions) for lack of anything better.

#NEW 1/2022 - found a polynomial parameterization of debye-waller factors.
#https://journals.iucr.org/a/issues/1999/05/00/sp0169/sp0169.pdf

#quoth KC Shambhu, I tried for Co, Fe, and Ge with reference to the values that we used for the high moment paper. Those values were taken at 295K. Using the parametrization, I get a 5% difference for Fe, 2 % for Ge, and 0.5% for Co. So, we can use this method to get B for Mn, which I get as 0.385 at 295K. [values he's citing are those noted below for those elements]


ScatteringFactor = [[0 for x in range(14)] for x in range(111)]  #roentgenium

    #Ti: B=0.55 https://www.publish.csiro.au/ph/pdf/ph880461

#0-4 are a_i, 5-9 are b_i, 10 is c in analytical expansion for fo 
#11, 12 = f', f". f' includes nuclear-Thompson and relativistic bits
#13 = Debye-Waller, for elemental xtal

ScatteringFactor[22][0] = 9.818524  #a1			#Ti
ScatteringFactor[22][5] = 8.001879   #b1
ScatteringFactor[22][1] = 1.522646   #a2
ScatteringFactor[22][6] = 0.029763   #b2
ScatteringFactor[22][2] = 1.703101   #a3
ScatteringFactor[22][7] = 39.885423  #b3
ScatteringFactor[22][3] = 1.768774  #a4
ScatteringFactor[22][8] = 120.1580  #b4
ScatteringFactor[22][4] = 7.082555   #a5
ScatteringFactor[22][9] = 0.532405   #b5
ScatteringFactor[22][10] = 0.102473 #c

if (XRAY=="Co"):
	ScatteringFactor[22][11] = -0.15370  #f'
	ScatteringFactor[22][12] = 2.3142*1j #f"
elif (XRAY=="Cu"):
	ScatteringFactor[22][11] = -0.13049
	ScatteringFactor[22][12] = 1.8070*1j
else: #default to Cu Ka if undefined
	ScatteringFactor[22][11] = -0.13049
	ScatteringFactor[22][12] = 1.8070*1j
ScatteringFactor[22][13] = 0.5173	#B 
#B from L.-M. Peng, G. Ren, S. L. Dudarev and M. J. Whelan @ 295K 

ScatteringFactor[26][0] = 12.311098  #a1			#Fe
ScatteringFactor[26][5] = 5.009415   #b1
ScatteringFactor[26][1] = 1.876623   #a2
ScatteringFactor[26][6] = 0.014461   #b2
ScatteringFactor[26][2] = 3.066177   #a3
ScatteringFactor[26][7] = 18.743041  #b3
ScatteringFactor[26][3] = 2.070451   #a4
ScatteringFactor[26][8] = 82.767874  #b4
ScatteringFactor[26][4] = 6.975185   #a5
ScatteringFactor[26][9] = 0.346506   #b5
ScatteringFactor[26][10] = -0.304931 #c

if (XRAY=="Co"):
	ScatteringFactor[26][11] = -3.3891  #f'
	ScatteringFactor[26][12] = 0.47507*1j #f"
elif (XRAY=="Cu"):
	ScatteringFactor[26][11] = -1.285
	ScatteringFactor[26][12] = 3.185*1j
else: #default to Cu Ka if undefined
	ScatteringFactor[26][11] = -1.285
	ScatteringFactor[26][12] = 3.185*1j
ScatteringFactor[26][13] = 0.3272	#B 
#B from L.-M. Peng, G. Ren, S. L. Dudarev and M. J. Whelan @ 295K bcc

ScatteringFactor[27][0] = 12.914510 #a1			#Co
ScatteringFactor[27][5] = 4.507138  #b1
ScatteringFactor[27][1] = 2.481908  #a2
ScatteringFactor[27][6] = 0.009126  #b2
ScatteringFactor[27][2] = 3.466894  #a3
ScatteringFactor[27][7] = 16.438130 #b3
ScatteringFactor[27][3] = 2.106351  #a4
ScatteringFactor[27][8] = 76.987317 #b4
ScatteringFactor[27][4] = 6.960892  #a5
ScatteringFactor[27][9] = 0.314418  #b5
ScatteringFactor[27][10] = -0.936572 #c
if (XRAY=="Co"):
	ScatteringFactor[27][11] = -2.0998  #f'
	ScatteringFactor[27][12] = 0.55705*1j #f"
elif (XRAY=="Cu"):
	ScatteringFactor[27][11] = -2.7647
	ScatteringFactor[27][12] = 3.6398*1j
else: #default to Cu Ka if undefined
	ScatteringFactor[27][11] = -2.7647
	ScatteringFactor[27][12] = 3.6398*1j
ScatteringFactor[27][13] = 0.307    #B	#https://onlinelibrary.wiley.com/iucr/doi/10.1107/S0108767399005176
#also a value of 0.39 at https://www.publish.csiro.au/ph/pdf/ph880461

ScatteringFactor[32][0] = 16.540614 #a1			#Ge
ScatteringFactor[32][5] = 2.866618  #b1
ScatteringFactor[32][1] = 1.567900  #a2
ScatteringFactor[32][6] = 0.012198  #b2
ScatteringFactor[32][2] = 3.727829  #a3
ScatteringFactor[32][7] = 13.432163 #b3
ScatteringFactor[32][3] = 3.345098  #a4
ScatteringFactor[32][8] = 58.866046 #b4
ScatteringFactor[32][4] = 6.785079  #a5
ScatteringFactor[32][9] = 0.210974  #b5
ScatteringFactor[32][10] = 0.018726 #c
if (XRAY=="Co"):
	ScatteringFactor[32][11] = -0.72563  #f'
	ScatteringFactor[32][12] = 1.1446*1j #f"
elif (XRAY=="Cu"):	
	ScatteringFactor[32][11] = -1.1475
	ScatteringFactor[32][12] = 0.88279*1j
else:	
	ScatteringFactor[32][11] = -1.1475
	ScatteringFactor[32][12] = 0.88279*1j
ScatteringFactor[32][13] = 0.6041		#B
#B from L.-M. Peng, G. Ren, S. L. Dudarev and M. J. Whelan @ 295K
	

# generate atomic scattering factor with the data in the matrices above
# we are not doing the thickness/absorption correction yet
# for films, you need to add the absorption/thickness correction later to each peak
# can turn off dispersion (f', f) and Debye-Waller corrections to compare with other 
# software or just see how much they matter. turning off dispersion also turns off
# the nuclear-Thompson and relativistic corrections to f' of course

def f(element,d):  			
	s = 1.0/(2.0*d)    #sin(theta)/lambda
	f=0
	for i in range(5):
		f += ScatteringFactor[element][i]*exp(-ScatteringFactor[element][i+5]*s*s)
		#sum a_i * exp(-b_i s^2)  gaussian approx to fo
	f += ScatteringFactor[element][10]  # + c 
	if (DISPERSION): #add dispersion f' and f''. note f'' is imaginary
		f += ScatteringFactor[element][11] + ScatteringFactor[element][12] 
	if (DEBYE_WALLER):
		f *= exp(-ScatteringFactor[element][13]*s*s)	
		#f_tot = (fo + f' + f'')*DW 
	return (f)		#note element [12] is imaginary, so return value is complex	

# general rules for a given space group on allowed hkl
			
def rules(h,k,l):	#general rules for allowed hkl  #are we handling permutable correctly?
	allowed=1
	if (space_group=="SG225" or space_group=="SG216"): #(225 and 216 have same rules)
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
	if (space_group=="SG224"):  
		if ((h==0) and ((k+l)%2==1)):
			allowed=0
		if ((k==0 and l==0) and (h%2==1)):
			allowed=0
		if (h==0 and k==0 and l==0):
			allowed=0			
	if (space_group=="SG194"):
		if ((h==k) and (l%2 == 1)):
			allowed=0
		if ((h==0 and k==0) and (l%2 == 1)):
			allowed=0
		if (h==0 and k==0 and l==0):
			allowed=0
	if (space_group=='SG139'):
		if ((h+k+l)%2==1):
			allowed = 0
		if ((h==0) and ((k+l)%2==1)):
			allowed = 0
		if ((l==0) and ((h+k)%2==1)):
			allowed = 0	
		if ((h==k) and (l%2==1)):
			allowed = 0
		if ((h==0 and k==0) and (l%2==1)):
			allowed = 0
		if ((k==0 and l==0) and (h%2==1)):
			allowed = 0
	return(allowed)

#calculate structure factor, including rules for specific sites in a space group

def F_hkl(site,h,k,l):		
	F=0
	if (space_group=="SG225" and (site[0]=='c8' or site[0]=='d24')): 
		if (h%2==0):
			for i in range (1,len(site)):
				F += exp(2*pi*1j*(site[i][0]*h+site[i][1]*k+site[i][2]*l))
	elif (space_group=="SG224" and (site[0]=='a2' or site[0]=='d6')): 
		if ((h+k+l)%2==0):
			for i in range (1,len(site)):
				F += exp(2*pi*1j*(site[i][0]*h+site[i][1]*k+site[i][2]*l))
	elif (space_group=="SG224" and (site[0]=='b4' or site[0]=='c4')): 
		if (((h+k)%2==0) and ((h+l)%2==0) and ((k+l)%2==0)):
			for i in range (1,len(site)):
				F += exp(2*pi*1j*(site[i][0]*h+site[i][1]*k+site[i][2]*l))	
	elif (space_group=="SG194" and (site[0]=='a2' or site[0]=='b2' or site[0]=='e4' or site[0]=='g6')):
		if (l%2==0):
			for i in range (1,len(site)):
				F += exp(2*pi*1j*(site[i][0]*h+site[i][1]*k+site[i][2]*l))
	elif (space_group=="SG194" and (site[0]=='c2' or site[0]=='d2' or site[0]=='f4')):
		if ((((l%2==0) or ((h-k)%3==1) or ((h-k)%3==2)))):
			for i in range (1,len(site)):
				F += exp(2*pi*1j*(site[i][0]*h+site[i][1]*k+site[i][2]*l))	
	elif (space_group=="SG139" and (site[0]=='c4' or site[0]=='d4')):
		if (l%2==0):	
			for i in range (1,len(site)):
				F += exp(2*pi*1j*(site[i][0]*h+site[i][1]*k+site[i][2]*l))	
	elif (space_group=="SG46"):
		if ((h!=0 and k!=0 and l!=0) and (h+k+l)%2==0):
			if (site[0]=='a4' and h%2!=0):
				F=0
			else: 
				for i in range (1,len(site)):
					F += exp(2*pi*1j*(site[i][0]*h+site[i][1]*k+site[i][2]*l))	
		if (h==0 and (k+l)%2==0):
			for i in range (1,len(site)):
				F += exp(2*pi*1j*(site[i][0]*h+site[i][1]*k+site[i][2]*l))		
		if (k==0 and (h+l)%2==0):
			if (site[0]=='a4' and h%2!=0):
				F=0
			else:
				for i in range (1,len(site)):
					F += exp(2*pi*1j*(site[i][0]*h+site[i][1]*k+site[i][2]*l))								
		if (l==0 and (h+k)%2==0):		
			if (site[0]=='a4' and h%2!=0):
				F=0
			else:
				for i in range (1,len(site)):
					F += exp(2*pi*1j*(site[i][0]*h+site[i][1]*k+site[i][2]*l))
		if (((h==0 and k==0) or (h==0 and l==0) or (k==0 and l==0))):
			if (h%2==0 and k%2==0 and l%2==0):
				for i in range (1,len(site)):
					F += exp(2*pi*1j*(site[i][0]*h+site[i][1]*k+site[i][2]*l))
	else: #any other site
		for i in range (1,len(site)):
				F += exp(2*pi*1j*(site[i][0]*h+site[i][1]*k+site[i][2]*l))
	return (F)

#find d spacing, bragg angle, and Lorentz-polarization factors

def d_hkl(h,k,l):								#d spacing for hkl
	if (space_group=="SG225" or space_group=="SG216" or space_group=="SG224"):
		tmp = (h**2+k**2+l**2)/(A*A)
	elif (space_group=="SG194"):
		tmp = (4.0/3.0)*(h*h+h*k+k*k)/(A*A) + l*l/(C*C)
	elif (space_group=="SG139"):
		tmp = (h**2+k**2)/(A*A) + l**2/(C*C)	
	elif (space_group=="SG46"):
		tmp = (h*h)/(A*A)+(k*k)/(B*B)+(l*l)/(C*C)
	return (1.0/sqrt(tmp))
	
def bragg(d,Lambda):   #just spits back 2theta given d and lambda 
	tmp = (Lambda/(2.0*d))
	if tmp<=1:
		angle = 2.0*degrees(arcsin(tmp)) 
	else: 
		angle = 0	#bad arcsin
	return(angle)   
	
def Lorentz_Pol(d,Lambda):      
	if ((Lambda/(2.0*d)) <= 1.0):
		if (SAMPLE_TYPE==POWDER):
			theta = arcsin(Lambda/(2.0*d)) #radians!
			LP = (1.0+(cos(2.0*theta))**2)/(sin(theta)*sin(2.0*theta))	
		elif (SAMPLE_TYPE==SINGLE_XTAL):
			theta = arcsin(Lambda/(2.0*d)) #radians!
			LP = (1.0+(cos(2.0*theta))**2)/(2.0*sin(2.0*theta))	
		else:
			theta = 0
			LP = 0			#bad specimen type
	else:
		theta = 0	#bad arcsin
		LP = 0
	return(LP)

#find all the peaks. no need for multiplicity factor since we brute force all combos
#yes, this is inelegant, but the code only takes a few seconds to run
			
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
				if (space_group=="SG194"):
					i = -(h+k) 	#fourth hexagonal index
				else:
					i = 0
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
					if (FILM):	#thickness correction factor
						G = 1.0-exp(-4.0*MU*THICKNESS*d/Lambda)
					else:
						G = 1.0 
					I = G*LP*(absolute(F_X1*f(X1[0],d)+F_Y1*f(Y1[0],d)+F_Z1*f(Z1[0],d)+F_X2*f(X2[0],d)+F_Y2*f(Y2[0],d)+F_Z2*f(Z2[0],d)))**2
					if ((two_theta<THETA_MAX) and (two_theta>THETA_MIN)):
						if (space_group=="SG194"): #if hex, output hkil
							pattern.append([two_theta,h,k,i,l,F_X1,F_X2,F_Y1,F_Y2,F_Z1,F_Z2,I,d])
						else: #if not hex, output hkl0
							pattern.append([two_theta,h,k,l,i,F_X1,F_X2,F_Y1,F_Y2,F_Z1,F_Z2,I,d])
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
		DataOut.writerow([space_group])
		DataOut.writerow(["a (A) lattice parameter"]+[A])
		if (space_group=="SG194" or space_group=="SG139"):
			DataOut.writerow(["c (A) lattice parameter"]+[C])
		DataOut.writerow(["Elements"])
		DataOut.writerow(["X1 X2 Y1 Y2 Z1 Z2 = "]+[(elements[str(X1[0])],elements[str(X2[0])],elements[str(Y1[0])],elements[str(Y2[0])],elements[str(Z1[0])],elements[str(Z2[0])])])		
		DataOut.writerow(["Sites"])
		DataOut.writerow(["X1 X2 Y1 Y2 Z1 Z2 = "]+[X1[1][0],X2[1][0],Y1[1][0],Y2[1][0],Z1[1][0],Z2[1][0]])
		DataOut.writerow(["Occupancy"])
		DataOut.writerow(["X1 X2 Y1 Y2 Z1 Z2 = "]+[X1[2],X2[2],Y1[2],Y2[2],Z1[2],Z2[2]])
	
		if (space_group=="SG194"):
			DataOut.writerow(['2T','h','k','i','l','Fx1','Fy1','Fz1','Fx2','Fy2','Fz2','I','d'])
		else:
			DataOut.writerow(['2T','h','k','l',' ','Fx1','Fy1','Fz1','Fx2','Fy2','Fz2','I','d'])
		DataOut.writerows(pattern)
		of.close()    	

	if (outputlistverbose):
		print("data for all allowed (hkl)")
		if (space_group=="SG194"):
			print("\n\n2Theta \t hkil \t X1 \t X2 \t Y1 \t Y2 \t Z1 \t Z2 \t I \t d (A)")
			for x in pattern:
				if x[11]!=0:
					print('{0:8.2f}\t ({1:1d},{2:1d},{3:1d},{4:1d}) \t {5:6.2f} \t {6:6.2f}  \t {7:8.2f} \t {8:6.2f} \t {9:6.2f}  \t {10:8.2f} \t {11}  \t {12:8.3f} '.format(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9],x[10],x[11],x[12]))
		else:
			print("\n\n2Theta \t hkl \t X1 \t X2 \t Y1 \t Y2 \t Z1 \t Z2 \t I \t d (A)")
			for x in pattern:
				if x[11]!=0:
					print('{0:8.2f}\t ({1:1d},{2:1d},{3:1d}) \t {4:6.2f} \t {5:6.2f}  \t {6:8.2f} \t {7:6.2f} \t {8:6.2f}  \t {9:8.2f} \t {10}  \t {11:8.3f} '.format(x[0], x[1], x[2], x[3], x[5], x[6], x[7], x[8], x[9],x[10],x[11],x[12]))

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
		print('\n2Theta \t I (normalized)')      #prints the results
		#TODO: add hkil list here - requires redoing dictionary
		#minor point as it is in the CSV file output

		for key,value in sorted(pattern_dict2.items()):
			if value!=0:
				print('{0:8.2f}\t{1:>10.6f}'.format(key,value))		

	if (outputfile):
		OutFile = "./output/"+elements[str(X1[0])]+X1[1][0]+elements[str(X2[0])]+X2[1][0]+elements[str(Y1[0])]+Y1[1][0]+elements[str(Y2[0])]+Y2[1][0]+elements[str(Z1[0])]+Z1[1][0]+elements[str(Z2[0])]+Z2[1][0]+"peak-list"+".csv"
		of = open(OutFile, 'wt')
		DataOut=csv.writer(of)
		DataOut.writerow(['2T','I (norm.)'])
		for key,value in sorted(pattern_dict2.items()):
			if value!=0:
				DataOut.writerow([key, value])
		of.close()    

	if (outputsites):	
		print("\n%s structure"%(space_group))
		print("Elements X1=%s X2=%s Y1=%s Y2=%s Z1=%s Z2=%s"%(elements[str(X1[0])],elements[str(X2[0])],elements[str(Y1[0])],elements[str(Y2[0])],elements[str(Z1[0])],elements[str(Z2[0])]))		
		print("Sites X1=%s X2=%s Y1=%s Y2=%s Z1=%s Z2=%s \nOccupancy X1=%s X2=%s Y1=%s Y2=%s Z1=%s Z2=%s \n "%(X1[1][0],X2[1][0],Y1[1][0],Y2[1][0],Z1[1][0],Z2[1][0],X1[2],X2[2],Y1[2],Y2[2],Z1[2],Z2[2]))	
		print("Max peak at %f\n"%maxtheta)	

	if (plot or plotfile):						#prepare data to plot
		x = [] 
		y = []
		for key,value in sorted(pattern_dict2.items()):
			x.append(key)
			if (plotsqrt):
				y.append(sqrt(value))
			else:
				y.append(value)
	
	if (plot):									#popup plot
		plt.bar(x,y,width=0.2,color='b',edgecolor='b')
		if (plotsqrt):
			plt.axis([THETA_MIN, THETA_MAX, 0, 10])
			plt.ylabel(r'$\sqrt{I}$'' (a.u.)')
		else:
			plt.axis([THETA_MIN, THETA_MAX, 0, 100])
			plt.ylabel('I (a.u.)')
		plt.xlabel(r'$2\theta\,(^\circ)$')
		plt.show()
	if (plotfile):								#just dump a PNG
		plt.bar(x,y,width=0.2,color='b',edgecolor='b')
		if (plotsqrt):
			plt.axis([THETA_MIN, THETA_MAX, 0, 10])
			plt.ylabel(r'$\sqrt{I}$'' (a.u.)')
		else:
			plt.axis([THETA_MIN, THETA_MAX, 0, 100])
			plt.ylabel('I (a.u.)')
		plt.xlabel(r'$2\theta\,(^\circ)$')
		foo = "./output/"+elements[str(X1[0])]+X1[1][0]+elements[str(X2[0])]+X2[1][0]+elements[str(Y1[0])]+Y1[1][0]+elements[str(Y2[0])]+Y2[1][0]+elements[str(Z1[0])]+Z1[1][0]+elements[str(Z2[0])]+Z2[1][0]+"graph"+".png"
		plt.savefig(foo, dpi = 600)
	return(maxtheta)
	
	
######## HERE IS A BLOCK THAT YOU EDIT ########	

#various output options
plot=1  					#pop-up plot of xrd pattern
plotfile=0					#save a png of the Pattern
plotsqrt=0					#use sqrt(I) for y axis in plots
outputfile=1				#output scattering factors & peak list to file
outputlist=1				#print a summary list of peaks to the tty
outputlistverbose=1			#print the ENTIRE list of peaks to the tty
outputsites=1				#print which elements are on which sites to tty

#discontinued functions, for now
#search_xyz=0				#search over x/y/z parameter 
#search_sites=0				#try switching site assignments [can take a while]

#Where are the elements? site X has [element,site,occupancy]
#if e.g., X1 and X2 elements share a site, take care that total occupancy isn't > 1

X1 = [Ti,Sites.c8,0.5]    		#2 atoms			#set as SG225 structure!
X2 = [Ti,Sites.c8,0.5] 			#2 atoms
Y1 = [Co,Sites.b4,0.5]			#1 atom
Y2 = [Co,Sites.b4,0.5]			#1 atom
Z1 = [Si,Sites.a4,0.5]			#1 atom
Z2 = [Si,Sites.a4,0.5]			#1 atom

# Sites.c8 = ['c8', (0,0,0),(0,0.5,0)]
# z=0
# Sites.a4 = ['a4', (0,0,0.2501)]
# 
# X1 = [Co,Sites.c8,1]   	
# X2 = [Co,Sites.a4,1] 		
# 
# Sites.b4 = ['b4', (0,0,0.5),(0,0.5,0.5)]
# 
# Y1 = [Ti,Sites.b4,0.5]
# Y2 = [Ti,Sites.b4,0.5]
# 
# Sites.c8 = ['c8', (0.25,0.25,0.25),(0.75,0.25,0.25)]
# Sites.b4 = ['b4', (0.25,0.9747,0.5055)]
# 
# Z1 = [Ga,Sites.c8,1]
# Z2 = [Ga,Sites.b4,1]    

#Sanity check: need X, Y, and Z total occupancies less than 1, not negative
#Partial occupancy is OK. But you know this, so we are just double checking.

#if ((X1[2]+X2[2]) > 1.0):
#	print("!!! Invalid X site occupancy, X1 + X2 > 1.")
#elif ((X1[2]<0) or (X2[2] < 0)):
#	print("!!! Invalid X site occupancy, X1 or X2 < 0.")
#elif ((Y1[2]+Y2[2]) > 1.0):
#	print("!!! Invalid Y site occupancy, Y1 + Y2 > 1.")
#elif ((Y1[2]<0) or (Y2[2] < 0)):
#	print("!!! Invalid Y site occupancy, Y1 or Y2 < 0.")	
#elif ((Z1[2]+Z2[2]) > 1.0):
#	print("!!! Invalid Z site occupancy, Z1 + Z2 > 1.")
#elif ((Z1[2]<0) or (Z2[2] < 0)):
#	print("!!! Invalid Z site occupancy, Z1 or Z2 < 0")	
#else:
Pattern(X1,X2,Y1,Y2,Z1,Z2,plot,outputfile,outputsites)	

X1tot = (len(X1[1])-1)*X1[2] 
X2tot = (len(X2[1])-1)*X2[2]
Y1tot = (len(Y1[1])-1)*Y1[2] 
Y2tot = (len(Y2[1])-1)*Y2[2]
Z1tot = (len(Z1[1])-1)*Z1[2] 
Z2tot = (len(Z2[1])-1)*Z2[2]

#Another sanity check - explicitly count the atoms up 

print("Composition used (X1 X2 Y1 Y2 Z1 Z2)")
print(elements[str(X1[0])],X1tot,elements[str(X2[0])],X2tot,elements[str(Y1[0])],Y1tot,elements[str(Y2[0])],Y2tot,elements[str(Z1[0])],Z1tot,elements[str(Z2[0])],Z2tot)


#run test patterns first so you know it is working

#if different atoms require different x, y, z for now you have to redefine sites as you go

#	Sites.c8 = ['c8', (x,y,z), (-x,y,z), (x+0.5,-y,z), (-x+0.5,y,z)]
#	Sites.b4 = ['b4', (0.25,y,z), (0.75,-y,z)]
#	Sites.a4 = ['a4', (0,0,z), (0.5,0,z)]


#following Jeitschko for FeTiSi
# x=0.5295
# y=0.1236
# Sites.c8 = ['c8', (x,y,z), (-x,y,z), (x+0.5,-y,z), (-x+0.5,y,z)]
# z=0
# Sites.a4 = ['a4', (0,0,z), (0.5,0,z)]
# 
# X1 = [Fe,Sites.c8,0.5]   
# X2 = [Fe,Sites.a4,0.5] 	
# 
# Sites.b4 = ['b4', (0.25,0.2207,0.0206), (0.25,0.4979,0.1677), (0.25,0.7996,0.0463), (0.75,-0.2207,0.0206), (0.75,-0.4979,0.1677), (0.27,-0.7996,0.0463)]
# 
# Y1 = [Ti,Sites.b4,1.0]
# Y2 = [Ti,Sites.b4,0]
# 
# x=0.506
# y=0.3325
# z=0.2452
# Sites.c8 = ['c8', (x,y,z), (-x,y,z), (x+0.5,-y,z), (-x+0.5,y,z)]
# y=0.0253
# z=0.2554
# Sites.b4 = ['b4', (0.25,y,z), (0.75,-y,z)]
# 
# Z1 = [Ge,Sites.c8,1.0]
# Z2 = [Ge,Sites.b4,1.0]    



#Example:vary Co content, but keep Ti+Co=2
# for i in range (0,11):
# 	x=i/10.0
# 	print"Co = 1+%f"%(x)
# 	X1 = [Fe,Sites.b4,1.0]    		#1 atom			#set as SG216 structure above!
# 	X2 = [Co,Sites.d4,1.0] 			#1 atom
# 	Y1 = [Co,Sites.c4,x]			#1 atom
# 	Y2 = [Ti,Sites.c4,1.0-x]		#1 atom
# 	Z1 = [Ge,Sites.a4,0.5]			#1 atom
# 	Z2 = [Ge,Sites.a4,0.5]			#1 atom
# 
# 	Pattern(X1,X2,Y1,Y2,Z1,Z2,plot,outputfile,outputsites)		
