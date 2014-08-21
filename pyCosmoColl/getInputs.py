'''
Code that collects the necessary inputs to the initial conditions code. Collected from screen prompts
'''

import sys
import os


def ICinputs():
	numerror = 0
	Mvir = None
	while not Mvir or Mvir < 0:
	        try: Mvir = float(raw_input("Virial Mass of halo at final redshift (float, e.g. 1e13):\n"))
	        except ValueError:
	                print 'invalid input, please try again'
	                numerror = 1
	        if numerror == 0 and Mvir <= 0: print "Invalid input, must be positive!"
	numerror = 0
	zstart = None
	while not zstart or zstart < 0:
	        try: zstart = float(raw_input("Initial redshift (float, e.g. 6):\n"))
	        except ValueError:
	                print 'invalid input, please try again'
	                numerror = 1
	        if numerror == 0 and zstart <= 0: print "Invalid input, must be positive!"
	zend = None
	numerror = 0
	while not zend or zend < 0:
	        try: zend = float(raw_input("Final redshift (float, e.g. 1):\n"))
	        except ValueError:
	                print 'invalid input, please try again'
	                numerror = 1
	        if numerror == 0 and zend < 0: print "Invalid input, must be positive!"
	        if zend == 0: break
	Ndark = None
	numerror = 0
	while not Ndark or Ndark < 0:
	        try: Ndark = float(raw_input("Number of dark matter particles (int, e.g. 10000):\n"))
	        except ValueError:
	                print 'invalid input, please try again'
	                numerror = 1
	        if numerror == 0 and Ndark <= 0: print "Invalid input, must be positive!"
	Ngas = None
	numerror = 0
	while not Ngas or Ngas < 0:
	        try: Ngas = float(raw_input("Number of gas particles (int, e.g. 10000):\n"))
	        except ValueError:
	                print 'invalid input, please try again'
	                numerror = 1
	        if numerror == 0 and Ngas < 0: print "Invalid input, must be positive!"
	        if Ngas == 0: break
	avir = None
	numerror = 0
	while not avir or avir < 0:
	        try: avir = float(raw_input("Fraction (1/x) of the initial sphere radius contining virialized mass by z = zend (float >1, e.g. 1.5):\n"))
	        except ValueError:
	                print 'invalid input, please try again'
	                numerror = 1
	        if numerror == 0 and avir < 1: print "Invalid input, must be > 1"
	Lambda = None
	numerror = 0
	while not Lambda or Lambda < 0 or Lambda > 1:
	        try: Lambda = float(raw_input("Angular Momentum Factor, Lambda [Lambda= L*E**0.5 / G*M**(5/2)], E ~ G*M**2/ R (float, e.g. 0.005) :\n"))
	        except ValueError:
	                print 'invalid input, please try again'
	                numerror = 1
	        if numerror == 0 and Lambda < 0: print "Invalid input, must be positive!"
	        if Lambda == 0: break
	KdivE = None
	numerror = 0
	while not KdivE or KdivE < 0:
	        try: KdivE = float(raw_input("Abs value of the Ratio of random tangential kinetic energy given to particles to their total energy. Nonzero value necessary (float, e.g. 0.001):\n"))
	        except ValueError:
	                print 'invalid input, please try again'
	                numerror = 1
	        if numerror == 0 and KdivE < 0: print "Invalid input, must be positive!"
	        if KdivE == 0: break
	fracH = None
	numerror = 0
	while not fracH or fracH < 0:
	        try: fracH = float(raw_input("maximum fraction of Hubble flow energy attained by random motion (float, e.g. 0.1):\n"))
	        except ValueError:
	                print 'invalid input, please try again'
	                numerror = 1
	        if numerror == 0 and fracH < 0: print "Invalid input, must be positive!"
	OverSamp = None
	numerror = 0
	while not OverSamp or OverSamp < 0:
	        try: OverSamp = float(raw_input("Factor by which dark matter is treated as 'oversampled'. Generally, this is equal to (Ndark/Ngas)^(1/3). This affects softening length calculation (fraction or Integer, e.g. 2, 3./2.):\n"))
	        except ValueError:
	                print 'invalid input, please try again'
       		        numerror = 1
        	if numerror == 0 and OverSamp < 0: print "Invalid input, must be positive!"
	temp = None
	if Ngas > 0:
        	numerror = 0
        	while not temp or temp < 0:
        	        try: temp = float(raw_input("Background Temperature,T_back, for adiabatic temperature profile T(r) = T_back * (1+delta)^(2/3) (float, e.g. 10000.):\n"))
        	        except ValueError:
        	                print 'invalid input, please try again'
        	                numerror = 1
        	        if numerror == 0 and temp < 0: print "Invalid input, must be positive!"	

	Name = str(raw_input("Name of file. Program will output files named YourName.bin and YourName.tbin (string):\n"))

	return Mvir, zstart, zend, Ndark, Ngas, avir,Lambda, KdivE,fracH,OverSamp,temp,Name
