import sys
from numpy import loadtxt

'''
Reads results collected with the reap.sh script into python.
Using IPython do:

%run -i getresults.py somename

where 'somename' is the prefix used when creating and collecting
the results, i.e. with dodirs.sh and reap.sh
'''

'''
params:  [Ek, m, symm, ricb, bci, bco, projection, forcing, \
			forcing_amplitude, forcing_frequency, magnetic, Em, Le2, N, lmax, toc1-tic, ncpus]
			
ken_dis: [KP, KT, internal_dis, rekin_dis, imkin_dis, repower, impower]
'''

if len(sys.argv) == 2:
	u = loadtxt(sys.argv[1]+'.flo')	# flow data
	p = loadtxt(sys.argv[1]+'.par')	# parameters
else:
	u = loadtxt('flow.dat')  		# flow data
	p = loadtxt('params.dat')  		# parameters
	
if len(u.shape)==1:
	u = u.reshape((-1,len(u)))
	p = p.reshape((-1,len(p)))
	

KP = u[:,0] 						# Poloidal kinetic energy
KT = u[:,1] 						# Toroidal kinetic energy
K  = KP + KT

ricb = p[:,3] 						# inner core radius
wf   = p[:,9] 						# forcing frequency
Ek   = p[:,0] 						# Ekman number

Dkin = u[:,3]*Ek 					# Kinetic energy dissipation
Dint = u[:,2]*Ek 					# Internal energy dissipation
rpow = u[:,5] 						# Input power from body forcing, real part
ipow = u[:,6] 						# Input power from body forcing, imaginary part
 
forcing = p[:,7]

bci = p[:,4]						# inner core boundary condition
bco = p[:,5]						# cmb boundary condition
amp = p[:,8]						# forcing amplitude

if shape(p)[1]>=15:
	N    = p[:,13]
	lmax = p[:,14]

if sum(forcing) == 0: 	# reads eigenvalue data
	w = loadtxt(sys.argv[1]+'.eig')
	if len(w.shape)==1:	
		w = w.reshape((-1,len(w)))

err1 = abs(-Dint/Dkin -1)

if shape(u)[1]>=10:
	vd1 = u[:,7]+u[:,8]
	vd2 = u[:,7]+u[:,9]
	
if shape(p)[1]>=17:
	t = p[:,15]
	ncpus = p[:,16]
   
magnetic = p[:,10]
if sum(magnetic) == np.shape(p)[0]: # reads magnetic data

	if len(sys.argv) == 2:
		b = loadtxt(sys.argv[1]+'.mag')
	else:
		b = loadtxt('magnetic.dat')
		
	if len(b.shape)==1:	
		b = b.reshape((-1,len(b)))
		
	M    = b[:,0] + b[:,1] 			# Total magnetic field energy
	Le2  = p[:,12] 					# Lehnert number squared
	Em   = p[:,11] 					# Magnetic Ekman number
	Dohm = (b[:,2]+b[:,3])*Le2*Em 	# Ohmic dissipation
	Le   = sqrt(Le2) 				# Lehnert number

	d = Dohm/Dint 					# Ohmic to viscous dissipation ratio
	
	if sum(Le) != 0:
		z = d*Em/Le**2

	if shape(b)[1]>=7 :
		od1 = b[:,4]+b[:,5]
		od2 = b[:,4]+b[:,6]
		d1 = od1/vd1				# dissip ratio in the bulk, without boundary layers
		d2 = od2/vd2				# a bit deeper in the bulk
		
		if sum(Le) != 0:
			z1 = d1*Em/Le**2
			z2 = d2*Em/Le**2
else:
	Dohm = 0
	
if sum(forcing) > 0:
	if sum(rpow) != 0:						# body forcing (input power should match total dissipation)
		err2 = abs( (rpow-(Dohm-Dkin))/rpow )
	else:									# boundary flow forcing (input power to be implemented)
		err2 = -1	
elif sum(forcing) == 0:						# eigenvalue problem (damping should match total dissipation)
	err2 = abs( 1+(Dohm-Dkin)/(w[:,0]*K) )

# total dissipation
D = Dint + Dohm

	

	
	
