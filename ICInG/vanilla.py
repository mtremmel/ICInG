import h5py
import struct
import numpy as np
import sys
import os
from pynbody.analysis import pkdgrav_cosmo as cosmo
import gc

def writetipsy(MGas, MDark, Ntot, nDark, nGas, posDark, posGas, velDark, velGas, temperature, hsml,gsoft,zstart,Name):
	'''
	Dump particle data to disk. Create a tipsy formatted file ".bin"
	Note: To use with ChaNGa or Gasoline, must convert to a tipsy binary using totipstd
	'''
        print "writing data to binary file. . ."
        file = []
        file.append(struct.pack('d i i i i i i', 1./(1+zstart), Ntot, 3, nGas, nDark, 0, 0))
        for i in np.arange(nGas):
            file.append(
                struct.pack('f f f f f f f f f f f f',
                            MGas[i], posGas[i,0], posGas[i,1], posGas[i,2],
                            velGas[i,0], velGas[i,1], velGas[i,2], 0., temperature[i],
                            hsml, 0., 0.,) #leave density, metals, and potential energy 0
            )
	for i in np.arange(nDark):
   		file.append(
   		     struct.pack('f f f f f f f f f',
   		                 MDark[i], posDark[i,0], posDark[i,1], posDark[i,2],
   		                 velDark[i,0], velDark[i,1], velDark[i,2], gsoft, 0.,) #leave potential energy 0
   		 )

	s = "".join(file)
	with open(Name+'.bin', 'w') as f:
	    f.write(s)
	
	print "converting binary..."
	os.system("totipstd <"+Name+".bin> "+Name+".tbin")
	os.system("rm "+Name+".bin")
	return


def delta(delta0,x,xf):   
	'''
	returns overdensity as a function of radius from center. Will be weird for x > xf so keep Rsphere <= xf. Taken from Evrard, 1988
	---inputs---
	delta0  = central overcensity
	x	= distance from center
	xf 	= size of sphere
	___output---
	overdinsity at radius x
	'''
        return (0.5*delta0)*(1.0+np.cos(x*np.pi/xf))

def perturb(x,delta0,xf):
	'''
	map constant density sphere to perturbed sphere. Evrard, 1988
	---inputs---
	x      = initial radial position
	delta0 = central over density
	xf     = size of sphere
	---output---
	perturbed position of particle with initial radial position x assuming isotropic sphere
	'''
        return 1 - delta(delta0,x,xf)/3.0
def Mr(x,delta0,xf,z):          
	'''
	Mass within radius x of a perturbed distribution from Evrard, 1988
	---inputs---
	x   = radial position
	delta0 = central over density
	xf = size of sphere
	z = redshift
	---output---
	M(< x)
	'''
        return (4.*np.pi*(1+z)**3) * ((x**3)/3. + (delta0/2.)*(x**3/3. + xf*((np.pi**2*x**2 - 2.*xf**2)*np.sin(np.pi*x/xf) + 2.*np.pi*x*xf*np.cos(np.pi*x/xf))/np.pi**3))
def Potential(r,delta0,rtot,z): 
	'''
	returns the gravitational potential given rho(r) = rho_crit(1+delta(r)), with delta(r) = 0.5*delta0*(1+cos(r*pi/R))
	---intputs---
	r = radial position
	delta0 = central over density
	rtot = radius of sphere
	z = redshift
	---output---
	phi(r)
	'''
        Minside = Mr(r,delta0,rtot,z)
        #solution to the integral of r*delta(r) from r to rtot
        Int_r_delta = 0.5*delta0*(0.5*(rtot**2-r**2) - (rtot*r/np.pi)*np.sin(np.pi*r/rtot) - (rtot**2/np.pi**2)*np.cos(np.pi*r/rtot) - (rtot/np.pi)**2)
        return - 1.0 * (Minside/r) - 4*np.pi*(1+z)**3*( (rtot**2 - r**2)/2. + Int_r_delta)


def IsoSphere(Ndark, Ngas,rsphere):
	'''
	creates an isotropic position for N particles using an input glass file, which is mirrored to create any number of particles
	-----inputs-----
	Ndark  = Number of DM particles
 	Ngas   = Number of gas particles
	rsphere = radius of the sphere you wish to create
 	-----outputs----
	posGas (position of Ngas particles)
	posDark (positions of Ndark particles)
	
	Nx3 dimensional arrays of the positions of all the N particles following the glass distribution.
	'''
	print "reading in glass data..."
	# Read glass data from file.  2**18 particles. Change line below for different init file
	_ROOT = os.path.abspath(os.path.dirname(__file__))
	gfname = 'data/snapshot_005.hdf5'
	glassfile = os.path.join(_ROOT, gfname)
	f = h5py.File(glassfile)
	coords = f['PartType1']['Coordinates'][:]
	x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
	del(coords)
	del(f)
	gc.collect()
	# Boxsize, from init glass file.
	bs = 1e5
	# See how many particles are in init file
	Nglass = len(x)
	# Mirror dataset multiple times  to generate enough particles
	print "mirroring data set. . ."
	while Nglass < 4*max(Ndark,Ngas):
	    x = np.concatenate([x, x, x-bs, x-bs, x, x, x-bs, x-bs]) + bs
	    y = np.concatenate([y, y-bs, y-bs, y, y, y-bs, y-bs, y]) + bs
	    z = np.concatenate([z, z, z, z, z-bs, z-bs, z-bs, z-bs]) + bs
	    bs *= 2
	    Nglass *= 8
	print "mirroring complete. Number of mirrors:", bs, Nglass, Ndark, Ngas, len(x)
	pos = np.transpose(np.array([x, y, z]))
	print "cleaning up..."
	del(x)
	del(y)
	del(z)
	gc.collect()
	print "sorting particles...."
	#reposition, sort by radius
	pos[:,0] -= bs/2
	pos[:,1] -= bs/2
	pos[:,2] -= bs/2
	rinit = np.sqrt(pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)
	radargs = np.argsort(rinit)
	posDark = np.array([])
	posGas = np.array([])
	rDark = 0
	rGas = 0
	print "setting final particle positions..."
	# Only accept the first n particles.  This won't work if npart ~ nglass.
	if Ndark > 0:
	        posDark = pos[radargs[0:Ndark],:]
	        rDark_init = rinit[radargs[Ndark-1]]
		posDark = posDark*(rsphere/rDark_init)
		rDark = np.sqrt(posDark[:,0]**2+posDark[:,1]**2+posDark[:,2]**2)
	gc.collect()
	if Ngas > 0:
	        posGas = pos[radargs[0:Ngas],:]
	        rGas_init = rinit[radargs[Ngas-1]]
		posGas = posGas*(rsphere/rGas_init)
		rGas = np.sqrt(posGas[:,0]**2+posGas[:,1]**2+posGas[:,2]**2)
	print "cleaning up again..."
	del(pos)
	del(rinit)
	del(radargs)
	gc.collect()
	return posDark, posGas, rDark, rGas

def genIC(CoMove=0,fmt = 'tipsy', Mvir=5e11,zstart=6,zend=0,Ndark=2**23,Ngas=0,avir=1.3,Lambda=0.05,TdivR=0.5,OverSamp=1,temp=0,Name='IsolatedCollapse'):
	'''
Code that creates a density perturbation IC under an Einstein-deSitter cosmology
Written by Michael Tremmel (2014)

Use of code:
(in ipython)
import ICInG
ICInG.vanilla.genIC()

requires tipsy_tools (see README file)

----Inputs----
CoMove  = 1 if user wants ICs in comoving coordinates (suggested: use 0)
fmt	= file format for output. ONLY SUPPORTED FOR TIPSY FORMAT RIGHT NOW
Mvir (in Solar masses) = final virial mass of the collapsed object required at the end of the simulation defined by zend
zstart  = starting redshift of the simulation. Suggestion: use 6
zend    = ending redshift of the simulation. Suggestion: use 0
Ndark   = number of dark matter particles desired (total gas+dark will be 2*Npart)
Ngas    = number of gas particles desired ("oversample" -> Ndark > Ngas)
avir    = the fraction of the initial sphirical radius that will be virialized at z = zend (rvir = Rsphere/avir). Actual radius is set by a combination of Mvir, a, zend, and zstart. Suggestion: use 1.3
lambda  = Factor that gives the total angular momentum of the system (lambda = L*E**0.5 / G*M**(5/2), E ~ gravitational binding energy ~ G*M**2/ R (not a bad approx for now...) Suggestion: use 0.05
TdivR     = Ratio of random tangential motion to peculiar radial velocity. Suggestion: use 0.5
OverSamp = factor by which dark matter has been "oversampled". e.g. 2 -> 2:1 -> 8 times the "normal" number of DM particles (i.e. normally corresponding to Ndm = 8*Ngas)
temp     = background temperature or the gas. Only matters if Ngas > 0. The code assumes that T(r) = temp * (1+delta(r))**(2/3), i.e. adiabatic contraction from a uniform temperature of temp.
Name    = Name desired for IC binary files (Name.bin, Name.tbin) 

----Output----
A tipsy standard format binary file (.tbin)
The Units will be cosmo50 units (see README file)

	'''
	#get necessary inputs
	
	print "The parameters you have chosen are:\n"
	print "Mvir: ", Mvir, "\n"
	print "zstart: ", zstart, "\n"
	print "zend: ", zend, "\n"
	print "Ndark: ", Ndark, "\n"
	print "Ngas: ", Ngas, "\n"
	print "avir: ", avir, "\n"
	print "Lambda: ", Lambda, "\n"
	print "TdivR: ", TdivR, "\n"
	print "OverSamp: ", OverSamp, "\n"
	print "temp: ", temp, "\n"
	print "Name: ", Name, "\n"

	Ntot = Ndark + Ngas
	#Units!
	dMsolUnit       = 1.84793e16
	dKpcUnit        = 50000.

	Mvir = Mvir / dMsolUnit
	#baryon fraction of mass in Universe (based on latest WMAP data)
	if Ngas > 0: fbar = 0.18181818
	if Ngas == 0: fbar = 0.0
	#Softening reference. Taken from DM only simualtion with Mdark = 1.68e5 Msun
	softRef = .173/dKpcUnit
	MpartRef = 1.68e5/dMsolUnit
	
	print "Calculating free parameters in density distribution"
	#set central overdensity peak using avir and requiring virialization by z = zend.
	overDenVir = 1.7*(1+zend)/(1+zstart)
	delta0 = (2.0*np.pi**3*overDenVir)  /  (3.0*(np.pi**2 - 2.0*avir**2) * avir*np.sin(np.pi/avir)   +   6.0*np.pi * avir**2 * np.cos(np.pi/avir)  +  np.pi**3)
	#set rsphere by requiring that Mvir must be inside rsphere/avir
	rsphere =((Mvir*3.0*avir**3*np.pi**2)  /  (2.0 * (1+zstart)**3 * (3.0*(np.pi**2 - 2.0*avir**2) * avir*delta0*np.sin(np.pi/avir)  +  6.0*np.pi * delta0*avir**2 * np.cos(np.pi/avir)  +  np.pi**3*(delta0 + 2))))**(1.0/3.0)
	
	
	posDark, posGas, rDark, rGas = IsoSphere(Ndark, Ngas,rsphere)

	print "perturbing distribution. . ."
	posDark_Perturb = np.array([])
	posGas_Perturb = np.array([])
	if Ndark > 0:
	        posDark_Perturb = np.zeros_like(posDark)
	        posDark_Perturb[:,0] = posDark[:,0] * perturb(rDark,delta0,rsphere)
	        posDark_Perturb[:,1] = posDark[:,1] * perturb(rDark,delta0,rsphere)
	        posDark_Perturb[:,2] = posDark[:,2] * perturb(rDark,delta0,rsphere)
	        rDark_Perturb = np.sqrt(posDark_Perturb[:,0]**2+posDark_Perturb[:,1]**2+posDark_Perturb[:,2]**2)
		del(posDark)
	if Ngas > 0:
	        posGas_Perturb = np.zeros_like(posGas)
	        posGas_Perturb[:,0] = posGas[:,0] * perturb(rGas,delta0,rsphere)
	        posGas_Perturb[:,1] = posGas[:,1] * perturb(rGas,delta0,rsphere)
	        posGas_Perturb[:,2] = posGas[:,2] * perturb(rGas,delta0,rsphere)
	        rGas_Perturb = np.sqrt(posGas_Perturb[:,0]**2+posGas_Perturb[:,1]**2+posGas_Perturb[:,2]**2)
		del(posGas)
	gc.collect()
	print "defining particle masses"
	#Divide up mass for Gas, Dark matter evenly for each particle of each type. Remember we only know for sure the mass of the material inside the future virialized regino (R < rsphere/avir)
	MDark = np.zeros(Ndark)
	if Ndark > 0: MDark[:] = (1-fbar)*Mvir / np.size(np.where(rDark_Perturb < rsphere/avir))
	MGas = np.zeros(Ngas)
	if Ngas > 0: MGas[:] = fbar*Mvir / np.size(np.where(rGas_Perturb < rsphere/avir))
	
	#assign perturbed velocities (based on Evrard 1988)
	H0 = np.sqrt(8.0*np.pi/3.0)  #due to the fact that our units are such that rho_crit = 1, G = 1
	H = H0*(1+zstart)**(3.0/2.0) #transform to starting redshift value
	
	velDark = np.array([])
	velGas = np.array([])	
	
	if Ndark > 0:           #Set the Hubble flow for each particle 
	        velDark = np.zeros((Ndark,3))
	        vInFall_Dark = (2.0/3.0)*delta(delta0,rDark,rsphere)*rDark*H
	        vHubDark =  rDark_Perturb*H
		del(rDark)
	if Ngas > 0:
	        velGas = np.zeros((Ngas,3))
	        vInFall_Gas = (2.0/3.0)*delta(delta0,rGas,rsphere)*rGas*H
	        vHubGas = rGas_Perturb*H
		del(rGas)
	gc.collect()
	if Ndark > 0:
		if np.size(np.where(vHubDark < vInFall_Dark)) > 0: print "uh oh, some DM particles have too large inward velocities!"
	if Ngas > 0:
		if np.size(np.where(vHubGas < vInFall_Gas)) > 0: print "uh oh, some GAS particles have too large inward velocities!"
	
	#Put in angular momentum -- using a lot of approximations
	if Lambda > 0:
	        print "putting in angular momentum, Lambda ", Lambda
	        Mtot = MDark.sum()
		EgGas = np.array([])
		EgDark = np.array([])
		EkhubDark = np.array([])
		EkhubGas = np.array([])
		if Ndark >0:
			EgDark = Potential(rDark_Perturb,delta0,rsphere,zstart)*MDark
			EkhubDark = 0.5*MDark*(vHubDark - vInFall_Dark)**2
		if Ngas > 0: 
			EgGas = Potential(rGas_Perturb,delta0,rsphere,zstart)*MGas
			EkhubGas = 0.5*MGas*(vHubGas - vInFall_Gas)**2
	        Eg = EgDark.sum()+EgGas.sum() #Binding energy of the (constant density) sphere (to order unity)
		Ek = EkhubGas.sum()+EkhubDark.sum()
		if Eg + Ek > 0: print "Warning: radial velocities too high!"
		E = np.abs(Eg + Ek)
	        L = Lambda*Mtot**(5./2.)/E**(0.5) #From Peebles book
	        I = (2./5.)*Mtot*rsphere**2  #approximate moment of inertia to be that of a uniform sphere... the perturbations should be small enough for this to be relatively accurate
	        w = L/I
	        print "Angular momentum: ", L, "angular speed: ", w
	        if Ngas > 0:
	                rxyGas = np.sqrt(posGas_Perturb[:,0]**2+posGas_Perturb[:,1]**2)
	                vphi_Gas = rxyGas * w
	                velGas[:,0] += vphi_Gas*(-1.0)*(posGas_Perturb[:,1]/rxyGas)
	                velGas[:,1] += vphi_Gas*(posGas_Perturb[:,0]/rxyGas)
	        if Ndark > 0:
	                rxyDark = np.sqrt(posDark_Perturb[:,0]**2+posDark_Perturb[:,1]**2)
	                vphi_Dark = rxyDark * w
	                velDark[:,0] += vphi_Dark*(-1.0)*(posDark_Perturb[:,1]/rxyDark)
	                velDark[:,1] += vphi_Dark*(posDark_Perturb[:,0]/rxyDark)

	#put in random motions for DM
	if TdivR > 0 and Ndark > 0:
	        print "putting in random tangential velocities with Krand/Kradial", TdivR
		vmag = TdivR *vInFall_Dark
	        Omega = np.random.random(Ndark) * 2.*np.pi  #random angle [0,2pi)
	        Theta = np.arccos(posDark_Perturb[:,2]/rDark_Perturb)
	        Phi = np.arctan(posDark_Perturb[:,1]/posDark_Perturb[:,0])
	        vTanRand = np.zeros_like(posDark_Perturb)
	        velDark[:,0] += vmag*(np.sin(Omega)*np.cos(Theta)*np.cos(Phi) - np.cos(Omega)*np.sin(Phi))
	        velDark[:,1] += vmag*(np.cos(Omega)*np.cos(Phi) + np.sin(Omega)*np.cos(Theta)*np.sin(Phi))
	        velDark[:,2] += vmag*(-1.0*np.sin(Omega)*np.sin(Theta))
	

	if Ndark > 0:
	        velDark_rad = vHubDark - vInFall_Dark #unlike Evrard 1988, we are working in PHYSICAL units. Also, take out the energy from tangential motion   
		if TdivR > 0: velDark_rad = np.sqrt(velDark_rad**2 - vmag**2)
       		velDark[:,0] += velDark_rad*np.sin(np.arccos(posDark_Perturb[:,2]/rDark_Perturb))*posDark_Perturb[:,0]/np.sqrt(posDark_Perturb[:,0]**2+posDark_Perturb[:,1]**2)
        	velDark[:,1] += velDark_rad*np.sin(np.arccos(posDark_Perturb[:,2]/rDark_Perturb))*posDark_Perturb[:,1]/np.sqrt(posDark_Perturb[:,0]**2+posDark_Perturb[:,1]**2)
        	velDark[:,2] += velDark_rad*(posDark_Perturb[:,2]/rDark_Perturb)
		del(vHubDark)
		del(vInFall_Dark)
		del(velDark_rad)
	if Ngas > 0:
	        velGas_rad =  vHubGas - vInFall_Gas
	        velGas[:,0] += velGas_rad*np.sin(np.arccos(posGas_Perturb[:,2]/rGas))*posGas_Perturb[:,0]/np.sqrt(posGas_Perturb[:,0]**2+posGas_Perturb[:,1]**2)
	        velGas[:,1] += velGas_rad*np.sin(np.arccos(posGas_Perturb[:,2]/rGas))*posGas_Perturb[:,1]/np.sqrt(posGas_Perturb[:,0]**2+posGas_Perturb[:,1]**2)
	        velGas[:,2] += velGas_rad*(posGask_Perturb[:,2]/rGas_Perturb)
		del(vHubGas)
		del(vInFall_Gas)
		del(velGas_rad)
	gc.collect()
	#set gas temperature
	temperature = 0
	if Ngas > 0: temperature = temp * (1.+delta(delta0,rGas_Perturb,rsphere))**(2./3.)
	#pretty cold, but not too cold. Cold enough to be << Tvir for any virialized halos we will be simulating here... hopefully.
	#I assume initially the gas has contracted adiabatically with the environment, so T ~ rho^(gamma-1) (gamma - 1 = 2/3 for neutral hydrogen)
	
	
	#find good values for gravitational softening. Compare DM particle mass to a reference case (values given above).
	#If DM is oversampled (Ndm > Ng) then pretend that it isn't (Ndm = Ng) when doing the calculation
	print "calculating softening length"
	partSep = (((4./3.)*np.pi*(rsphere/avir)**3)/np.size(np.where(rDark_Perturb < rsphere/avir)))**(1./3.)
	if OverSamp <= 0:
	        print "WARNING: 0 or negative value for OverSamp. Assuming you mean Oversampled is 1. Must be > 0..."
	        OverSamp = 1
	gsoft = (MDark[0]/MpartRef)**(1./3.)*OverSamp*softRef
	if gsoft  > partSep:
	        gsoft = partSep
	
	#just set SPH smoothing to gravitational softening at beginning since it doesn't really matter.
	hsml = gsoft

	t0 = 2./(3.*H0)
        t = t0/(1+zstart)**(3./2)
        total_time = t0 - t	
	# Change to CoMoving if prompted
	if CoMove==1:
		aFac = (1+zstart)**-1
		if Ndark>0:
			velDark = (velDark - H*posDark_Perturb)/aFac
			posDark_Perturb = posDark_Perturb/aFac
		if Ngas>0:
                        velGas = (velGas - H*posGas_Perturb)/aFac
                        posGas_Perturb = posGas_Perturb/aFac
		
	#Dump to a tipsy format file

	if fmt == 'tipsy': writetipsy(MGas, MDark, Ntot, Ndark, Ngas, posDark_Perturb, posGas_Perturb, velDark, velGas, temperature, hsml, gsoft,zstart,Name)
	else: 
		print "error, filetype not regonized. Attempting to write tipsy format anyway"
		writetipsy(MGas, MDark, Ntot, Ndark, Ngas, posDark_Perturb, posGas_Perturb, velDark, velGas, temperature, hsml, gsoft,zstart,Name)
	#print out important things for user reference
	print "tbin file created in standard tipsy format!"
	if Ngas > 0: print "Mass per particle: ", MDark[0], "(DM)", MGas[0], "(gas)"
	else: print "Mass per particle: ", MDark[0], "(DM)"
	print "SPH soft ", hsml, " grav soft ", gsoft
	if Ngas > 0: print "Sphere size ", rsphere, " Sphere Mass ", np.sum(MDark)+np.sum(MGas)
	else: print "Sphere size (physical units)", rsphere, " Sphere Mass ", np.sum(MDark)
	print "delta0", delta0
	print "runtime (in sim units) for Z = 0", total_time
	if Lambda > 0: print "angular speed ", w	

	return






