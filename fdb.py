'''
   FDB.py
   Fast generation of FDB approach to compute E-field dependent chemical barriers

		Author:      PhD candidate Pau Besalú-Sala
		Affiliation: Institut de Química Computacional i Catálisi
			     Universitat de Girona
	        Year: 	     2021
	        Version:     FDB1.1

		Usage:	     See manual (FDB_py_manual.txt) 
'''

#charging modules
import numpy as np
import os

#Initialization
E0       = np.zeros(2)
mu       = np.zeros((2,3))
alpha    = np.zeros((2,3,3))
alphanr  = np.zeros((2,3,3))
beta     = np.zeros((2,3,3,3))
kcal	 = 6.2751e2

dmu      = np.zeros(3)
dalpha   = np.zeros((3,3))
dalphanr = np.zeros((3,3))
dbeta	 = np.zeros((3,3,3))

#Checking if all the files exist
check1 = os.path.isfile('fdb.inp')
check2 = os.path.isfile('reactants.dat')
check3 = os.path.isfile('ts.dat')
if check1 == False or check2 == False or check3 == False :
   print()
   if check1 == False :
     print(" Missing fdb.inp file in execution directory")
   if check2 == False :
     print(" Missing reactants.dat file in execution directory")
   if check3 == False :
     print(" Missing ts.dat file in execution directory")

   print(" Stop execution of python script.")
   print()
   quit()


#Opening the fdb.inp file
file1 = open('fdb.inp', 'r')
d1_aprox	= int(file1.readline())   #Approximation used for 1D

### Start of 1D FDB data generation ###
if d1_aprox == 1 or d1_aprox == 2 or d1_aprox == 3 or d1_aprox == 4:
   print()
   print(" 1D data generation ")
   print(" ------------------ ")
   print()

   # Open the file from which the properties are read
   file2 = open('reactants.dat', 'r' )
   file3 = open('ts.dat','r')
   E0[0] = float(file2.readline()) #The free field energy for R
   E0[1] = float(file3.readline()) #The free field energy for TS
   dE=E0[1]-E0[0]
   print(" Free-field Electronic energy barrier is %.5f kcal/mol " %(dE*kcal) )

   if d1_aprox == 1:
     print(" Aproximation used is only dipole")
     print()
     #Reactants data
     mu[0,0] = float(file2.readline()) #The dipole moment x component 
     mu[0,1] = float(file2.readline()) #The dipole moment y component
     mu[0,2] = float(file2.readline()) #The dipole moment z component
     #TS data
     mu[1,0] = float(file3.readline()) #The dipole moment x component 
     mu[1,1] = float(file3.readline()) #The dipole moment y component
     mu[1,2] = float(file3.readline()) #The dipole moment z component

     print(" Dipole moment for R and TS (au) are ")
     print(" x= %.5f %.5f " %(mu[0,0],mu[1,0]))
     print(" y= %.5f %.5f " %(mu[0,1],mu[1,1]))
     print(" z= %.5f %.5f " %(mu[0,2],mu[1,2]))
     print()

   elif d1_aprox == 2:
     print(" Aproximation used is dipole + polarizability(el)")
     print()
     #Reactants data
     mu[0,0] = float(file2.readline()) #The dipole moment x component
     mu[0,1] = float(file2.readline()) #The dipole moment y component
     mu[0,2] = float(file2.readline()) #The dipole moment z component
     alpha[0,0,0] = float(file2.readline()) # Polarizability xx component
     alpha[0,1,0] = float(file2.readline()) # Polarizability yx component
     alpha[0,1,1] = float(file2.readline()) # Polarizability yy component
     alpha[0,2,0] = float(file2.readline()) # Polarizability zx component
     alpha[0,2,1] = float(file2.readline()) # Polarizability zy component
     alpha[0,2,2] = float(file2.readline()) # Polarizability zz component

     #TS data
     mu[1,0] = float(file3.readline()) #The dipole moment x component
     mu[1,1] = float(file3.readline()) #The dipole moment y component
     mu[1,2] = float(file3.readline()) #The dipole moment z component
     alpha[1,0,0] = float(file3.readline()) # Polarizability xx component
     alpha[1,1,0] = float(file3.readline()) # Polarizability yx component
     alpha[1,1,1] = float(file3.readline()) # Polarizability yy component
     alpha[1,2,0] = float(file3.readline()) # Polarizability zx component
     alpha[1,2,1] = float(file3.readline()) # Polarizability zy component
     alpha[1,2,2] = float(file3.readline()) # Polarizability zz component

     print(" Dipole moment for R and TS (au) is ")
     print(" x= %.5f %.5f " %(mu[0,0],mu[1,0]))
     print(" y= %.5f %.5f " %(mu[0,1],mu[1,1]))
     print(" z= %.5f %.5f " %(mu[0,2],mu[1,2]))
     print()
     print(" Polarizability for R and TS (au) is ")
     print(" xx= %.5f %.5f" % (alpha[0,0,0],alpha[1,0,0]))
     print(" yx= %.5f %.5f" % (alpha[0,1,0],alpha[1,1,0]))
     print(" yy= %.5f %.5f" % (alpha[0,1,1],alpha[1,1,1]))
     print(" zx= %.5f %.5f" % (alpha[0,2,0],alpha[1,2,0]))
     print(" zy= %.5f %.5f" % (alpha[0,2,1],alpha[1,2,1]))
     print(" zz= %.5f %.5f" % (alpha[0,2,2],alpha[1,2,2]))
     print(" Polarizability matrix is considered symmetric.")
     print()

   elif d1_aprox == 3:
     print(" Aproximation used is dipole + polarizability(el) + 1st hyperpolarizability")
     print(" WARNING: Beta is changed of sign because G16 has an ill-definition.")
     print(" If you do not work with G16 please enter the input with changed sign for beta. ")
     print()
     #Reactants data
     mu[0,0] = float(file2.readline()) #The dipole moment x component
     mu[0,1] = float(file2.readline()) #The dipole moment y component
     mu[0,2] = float(file2.readline()) #The dipole moment z component
     alpha[0,0,0] = float(file2.readline())  # Polarizability xx component
     alpha[0,1,0] = float(file2.readline())  # Polarizability yx component
     alpha[0,1,1] = float(file2.readline())  # Polarizability yy component
     alpha[0,2,0] = float(file2.readline())  # Polarizability zx component
     alpha[0,2,1] = float(file2.readline())  # Polarizability zy component
     alpha[0,2,2] = float(file2.readline())  # Polarizability zz component
     beta[0,0,0,0] = -float(file2.readline()) # Beta xxx component 
     beta[0,0,0,1] = -float(file2.readline()) # Beta xxy component
     beta[0,1,0,1] = -float(file2.readline()) # Beta yxy component
     beta[0,1,1,1] = -float(file2.readline()) # Beta yyy component
     beta[0,0,0,2] = -float(file2.readline()) # Beta xxz component
     beta[0,1,0,2] = -float(file2.readline()) # Beta yxz component
     beta[0,1,1,2] = -float(file2.readline()) # Beta yyz component
     beta[0,2,0,2] = -float(file2.readline()) # Beta zxz component
     beta[0,2,1,2] = -float(file2.readline()) # Beta zyz component
     beta[0,2,2,2] = -float(file2.readline()) # Beta zzz component

     #TS data
     mu[1,0] = float(file3.readline()) #The dipole moment x component
     mu[1,1] = float(file3.readline()) #The dipole moment y component
     mu[1,2] = float(file3.readline()) #The dipole moment z component
     alpha[1,0,0] = float(file3.readline())  # Polarizability xx component
     alpha[1,1,0] = float(file3.readline())  # Polarizability yx component
     alpha[1,1,1] = float(file3.readline())  # Polarizability yy component
     alpha[1,2,0] = float(file3.readline())  # Polarizability zx component
     alpha[1,2,1] = float(file3.readline())  # Polarizability zy component
     alpha[1,2,2] = float(file3.readline())  # Polarizability zz component
     beta[1,0,0,0] = -float(file3.readline()) # Beta xxx component
     beta[1,0,0,1] = -float(file3.readline()) # Beta xxy component
     beta[1,1,0,1] = -float(file3.readline()) # Beta yxy component
     beta[1,1,1,1] = -float(file3.readline()) # Beta yyy component
     beta[1,0,0,2] = -float(file3.readline()) # Beta xxz component
     beta[1,1,0,2] = -float(file3.readline()) # Beta yxz component
     beta[1,1,1,2] = -float(file3.readline()) # Beta yyz component
     beta[1,2,0,2] = -float(file3.readline()) # Beta zxz component
     beta[1,2,1,2] = -float(file3.readline()) # Beta zyz component
     beta[1,2,2,2] = -float(file3.readline()) # Beta zzz component

     print(" Dipole moment for R and TS (au) is ")
     print(" x= %.5f %.5f " %(mu[0,0],mu[1,0]))
     print(" y= %.5f %.5f " %(mu[0,1],mu[1,1]))
     print(" z= %.5f %.5f " %(mu[0,2],mu[1,2]))
     print()
     print(" Polarizability for R and TS (au) is ")
     print(" xx= %.5f %.5f" % (alpha[0,0,0],alpha[1,0,0]))
     print(" yx= %.5f %.5f" % (alpha[0,1,0],alpha[1,1,0]))
     print(" yy= %.5f %.5f" % (alpha[0,1,1],alpha[1,1,1]))
     print(" zx= %.5f %.5f" % (alpha[0,2,0],alpha[1,2,0]))
     print(" zy= %.5f %.5f" % (alpha[0,2,1],alpha[1,2,1]))
     print(" zz= %.5f %.5f" % (alpha[0,2,2],alpha[1,2,2]))
     print(" Polarizability matrix is considered symmetric.")
     print()
     print(" 1st Hyperpolarizability for R and TS (au) is ")
     print(" xxx= %.5f %.5f " %(beta[0,0,0,0],beta[1,0,0,0]))
     print(" xxy= %.5f %.5f " %(beta[0,0,0,1],beta[1,0,0,1]))
     print(" yxy= %.5f %.5f " %(beta[0,1,0,1],beta[1,1,0,1]))
     print(" yyy= %.5f %.5f " %(beta[0,1,1,1],beta[1,1,1,1]))
     print(" xxz= %.5f %.5f " %(beta[0,0,0,2],beta[1,0,0,2]))
     print(" yxz= %.5f %.5f " %(beta[0,1,0,2],beta[1,1,0,2]))
     print(" yyz= %.5f %.5f " %(beta[0,1,1,2],beta[1,1,1,2]))
     print(" zxz= %.5f %.5f " %(beta[0,2,0,2],beta[1,2,0,2]))
     print(" zyz= %.5f %.5f " %(beta[0,2,1,2],beta[1,2,1,2]))
     print(" zzz= %.5f %.5f " %(beta[0,2,2,2],beta[1,2,2,2]))
     print(" Hyperpolarizability tensor is considered symmetric.")
     print()

   else: # therefore d1_aprox == 4
     print(" Aproximation used is dipole + polarizability(el+nr) + 1st hyperpolarizability")
     print(" WARNING: Beta is changed of sign because G16 has an ill-definition.")
     print(" If you do not work with G16 please enter the input with changed sign for beta. ")
     print()
     #Reactants data
     mu[0,0] = float(file2.readline()) #The dipole moment x component
     mu[0,1] = float(file2.readline()) #The dipole moment y component
     mu[0,2] = float(file2.readline()) #The dipole moment z component
     alpha[0,0,0] = float(file2.readline())  # Polarizability xx component
     alpha[0,1,0] = float(file2.readline())  # Polarizability yx component
     alpha[0,1,1] = float(file2.readline())  # Polarizability yy component
     alpha[0,2,0] = float(file2.readline())  # Polarizability zx component
     alpha[0,2,1] = float(file2.readline())  # Polarizability zy component
     alpha[0,2,2] = float(file2.readline())  # Polarizability zz component
     alphanr[0,0,0] = float(file2.readline())  # Polarizability xx component
     alphanr[0,1,0] = float(file2.readline())  # Polarizability yx component
     alphanr[0,1,1] = float(file2.readline())  # Polarizability yy component
     alphanr[0,2,0] = float(file2.readline())  # Polarizability zx component
     alphanr[0,2,1] = float(file2.readline())  # Polarizability zy component
     alphanr[0,2,2] = float(file2.readline())  # Polarizability zz component
     beta[0,0,0,0] = -float(file2.readline()) # Beta xxx component
     beta[0,0,0,1] = -float(file2.readline()) # Beta xxy component
     beta[0,1,0,1] = -float(file2.readline()) # Beta yxy component
     beta[0,1,1,1] = -float(file2.readline()) # Beta yyy component
     beta[0,0,0,2] = -float(file2.readline()) # Beta xxz component
     beta[0,1,0,2] = -float(file2.readline()) # Beta yxz component
     beta[0,1,1,2] = -float(file2.readline()) # Beta yyz component
     beta[0,2,0,2] = -float(file2.readline()) # Beta zxz component
     beta[0,2,1,2] = -float(file2.readline()) # Beta zyz component
     beta[0,2,2,2] = -float(file2.readline()) # Beta zzz component

     #TS data
     mu[1,0] = float(file3.readline()) #The dipole moment x component
     mu[1,1] = float(file3.readline()) #The dipole moment y component
     mu[1,2] = float(file3.readline()) #The dipole moment z component
     alpha[1,0,0] = float(file3.readline())  # Polarizability xx component
     alpha[1,1,0] = float(file3.readline())  # Polarizability yx component
     alpha[1,1,1] = float(file3.readline())  # Polarizability yy component
     alpha[1,2,0] = float(file3.readline())  # Polarizability zx component
     alpha[1,2,1] = float(file3.readline())  # Polarizability zy component
     alpha[1,2,2] = float(file3.readline())  # Polarizability zz component
     alphanr[1,0,0] = float(file3.readline())  # Polarizability xx component
     alphanr[1,1,0] = float(file3.readline())  # Polarizability yx component
     alphanr[1,1,1] = float(file3.readline())  # Polarizability yy component
     alphanr[1,2,0] = float(file3.readline())  # Polarizability zx component
     alphanr[1,2,1] = float(file3.readline())  # Polarizability zy component
     alphanr[1,2,2] = float(file3.readline())  # Polarizability zz component
     beta[1,0,0,0] = -float(file3.readline()) # Beta xxx component
     beta[1,0,0,1] = -float(file3.readline()) # Beta xxy component
     beta[1,1,0,1] = -float(file3.readline()) # Beta yxy component
     beta[1,1,1,1] = -float(file3.readline()) # Beta yyy component
     beta[1,0,0,2] = -float(file3.readline()) # Beta xxz component
     beta[1,1,0,2] = -float(file3.readline()) # Beta yxz component
     beta[1,1,1,2] = -float(file3.readline()) # Beta yyz component
     beta[1,2,0,2] = -float(file3.readline()) # Beta zxz component
     beta[1,2,1,2] = -float(file3.readline()) # Beta zyz component
     beta[1,2,2,2] = -float(file3.readline()) # Beta zzz component

     print(" Dipole moment for R and TS (au) is ")
     print(" x= %.5f %.5f " %(mu[0,0],mu[1,0]))
     print(" y= %.5f %.5f " %(mu[0,1],mu[1,1]))
     print(" z= %.5f %.5f " %(mu[0,2],mu[1,2]))
     print()
     print(" Polarizability for R and TS (au) is ")
     print(" xx= %.5f %.5f" % (alpha[0,0,0],alpha[1,0,0]))
     print(" yx= %.5f %.5f" % (alpha[0,1,0],alpha[1,1,0]))
     print(" yy= %.5f %.5f" % (alpha[0,1,1],alpha[1,1,1]))
     print(" zx= %.5f %.5f" % (alpha[0,2,0],alpha[1,2,0]))
     print(" zy= %.5f %.5f" % (alpha[0,2,1],alpha[1,2,1]))
     print(" zz= %.5f %.5f" % (alpha[0,2,2],alpha[1,2,2]))
     print()
     print(" Polarizability(nr) for R and TS (au) is ")
     print(" xx= %.5f %.5f" % (alphanr[0,0,0],alphanr[1,0,0]))
     print(" yx= %.5f %.5f" % (alphanr[0,1,0],alphanr[1,1,0]))
     print(" yy= %.5f %.5f" % (alphanr[0,1,1],alphanr[1,1,1]))
     print(" zx= %.5f %.5f" % (alphanr[0,2,0],alphanr[1,2,0]))
     print(" zy= %.5f %.5f" % (alphanr[0,2,1],alphanr[1,2,1]))
     print(" zz= %.5f %.5f" % (alphanr[0,2,2],alphanr[1,2,2]))
     print(" Polarizability matrixies are considered symmetric.")
     print()
     print(" 1st Hyperpolarizability for R and TS (au) is ")
     print(" xxx= %.5f %.5f " %(beta[0,0,0,0],beta[1,0,0,0]))
     print(" xxy= %.5f %.5f " %(beta[0,0,0,1],beta[1,0,0,1]))
     print(" yxy= %.5f %.5f " %(beta[0,1,0,1],beta[1,1,0,1]))
     print(" yyy= %.5f %.5f " %(beta[0,1,1,1],beta[1,1,1,1]))
     print(" xxz= %.5f %.5f " %(beta[0,0,0,2],beta[1,0,0,2]))
     print(" yxz= %.5f %.5f " %(beta[0,1,0,2],beta[1,1,0,2]))
     print(" yyz= %.5f %.5f " %(beta[0,1,1,2],beta[1,1,1,2]))
     print(" zxz= %.5f %.5f " %(beta[0,2,0,2],beta[1,2,0,2]))
     print(" zyz= %.5f %.5f " %(beta[0,2,1,2],beta[1,2,1,2]))
     print(" zzz= %.5f %.5f " %(beta[0,2,2,2],beta[1,2,2,2]))
     print(" Hyperpolarizability tensor is considered symmetric.")
     print()

   file2.close()
   file3.close()

   #Computing diferences on properties
   dmu      = mu[1,:]-mu[0,:]
   dalpha   = alpha[1,:,:]-alpha[0,:,:]
   dalphanr = alphanr[1,:,:]-alphanr[0,:,:]
   dbeta    = beta[1,:,:,:]-beta[0,:,:,:]

   numF = int(file1.readline())   #Number of EF desired to study
   F    = np.zeros(numF)
   FDBx = np.zeros(numF)
   FDBy = np.zeros(numF)
   FDBz = np.zeros(numF)
   for i in range(numF):
      F[i] = float(file1.readline()) #Each field to study (in a.u.)
   print(" Fields analyzed are ", F)
   print()

   # Calculation of each FDB
   if d1_aprox == 1 :
     for i in range(numF):
        FDBx[i]=(dE-dmu[0]*F[i])*kcal
        FDBy[i]=(dE-dmu[1]*F[i])*kcal
        FDBz[i]=(dE-dmu[2]*F[i])*kcal

   elif d1_aprox == 2 :
     for i in range(numF):
        FDBx[i]=(dE-dmu[0]*F[i]-0.5e0*dalpha[0,0]*F[i]**2)*kcal
        FDBy[i]=(dE-dmu[1]*F[i]-0.5e0*dalpha[1,1]*F[i]**2)*kcal
        FDBz[i]=(dE-dmu[2]*F[i]-0.5e0*dalpha[2,2]*F[i]**2)*kcal

   elif d1_aprox == 3 :
     for i in range(numF):
        FDBx[i]=(dE-dmu[0]*F[i]-0.5e0*dalpha[0,0]*F[i]**2-1.0e0/6.0e0*dbeta[0,0,0]*F[i]**3)*kcal
        FDBy[i]=(dE-dmu[1]*F[i]-0.5e0*dalpha[1,1]*F[i]**2-1.0e0/6.0e0*dbeta[1,1,1]*F[i]**3)*kcal
        FDBz[i]=(dE-dmu[2]*F[i]-0.5e0*dalpha[2,2]*F[i]**2-1.0e0/6.0e0*dbeta[2,2,2]*F[i]**3)*kcal

   else : # Therefore, d1_aprox=4
     for i in range(numF):
        FDBx[i]=(dE-dmu[0]*F[i]-0.5e0*(dalpha[0,0]+dalphanr[0,0])*F[i]**2-1.0e0/6.0e0*dbeta[0,0,0]*F[i]**3)*kcal
        FDBy[i]=(dE-dmu[1]*F[i]-0.5e0*(dalpha[1,1]+dalphanr[1,1])*F[i]**2-1.0e0/6.0e0*dbeta[1,1,1]*F[i]**3)*kcal
        FDBz[i]=(dE-dmu[2]*F[i]-0.5e0*(dalpha[2,2]+dalphanr[2,2])*F[i]**2-1.0e0/6.0e0*dbeta[2,2,2]*F[i]**3)*kcal

   ## Echoing the results
   print(" FINAL RESULTS")
   print(" -------------")
   print()
   print("  Field     FDBx         FDBy         FDBz ")
   print(" ------------------------------------------")
   for i in range(numF) :
     print("  %.5f  %.5f %.5f %.5f" %(F[i],FDBx[i],FDBy[i],FDBz[i]))
   print()

else:
   print()
   print(" Bad input for FDB approxmation in 1D data generation. Stop. ")
   quit()

### End of 1D FDB data generation ###


### Start of 2D FDB data generation ### 
d2_aprox = int(file1.readline()) # Approximation used for 2D

if d2_aprox == 0:
   print()
   print(" No 2D data generation. End of FDB.py run.")
   print()
   quit()

elif d2_aprox == 1 or d2_aprox == 2 or d2_aprox == 3 or d2_aprox == 4: 
   print()
   print(" 2D data generation ")
   print(" ------------------ ")
   print()
   if d2_aprox == 1:
     print(" Aproximation used is only dipole")
     print()
   elif d2_aprox == 2:
     print(" Aproximation used is dipole + polarizability(el)")
     print()
   elif d2_aprox == 3:
     print(" Aproximation used is dipole + polarizability(el) + 1st hyperpolarizability")
     print()
   else: #d2_aprox == 4
     print(" Aproximation used is dipole + polarizability(el+nr) + 1st hyperpolarizability")
     print()

   axis = int(file1.readline())    # Axis to be used, either xy, xz or yz
   fmin = float(file1.readline())  # Minimum field-strenght used (au)
   step = float(file1.readline())    # Density of the grid to be used (stepsize)
   fpoints = int(file1.readline()) # number of fields to analyze

   #Initializing
   #At this moment, fdb.py cannot handle different spans for the two axis.
   barrier_2D = np.zeros((fpoints,fpoints))
   nofield_2D = np.ones((fpoints,fpoints))
   nofield_2D = dE*nofield_2D
   f_print = np.zeros(fpoints) #Point to analyze (2D)
   mu_contrib      = np.zeros((fpoints,fpoints))
   alpha_contrib   = np.zeros((fpoints,fpoints))
   alphanr_contrib = np.zeros((fpoints,fpoints))
   beta_contrib    = np.zeros((fpoints,fpoints))
 
   for i,val in enumerate (f_print): 
      f_print[i] = fmin+step*i

   #Displaying some basic info
   print()
   print(" Minimum field-strenght (a.u.) is ", fmin)
   print(" Stepsize is ", step)
   print(" Number of points to analize is ", fpoints)
   print(" Therefore, maximum field-strenght (a.u.) is ", f_print[fpoints-1])
   print()
   
   #Do the actual calculation
   if axis == 1: #XY axis --> 0 and 1 labels 
      print(" 2D printing, xy axis")
      print()
      
      fz = 0.0e0
      for i,fx in enumerate (f_print):
       for j,fy in enumerate (f_print):
        mu_contrib[i,j]      = -1.0e0*(dmu[0]*fx+dmu[1]*fy+dmu[2]*fz)

        p = dalpha[0,0]*fx**2+dalpha[1,1]*fy**2+dalpha[2,2]*fz**2
        p = p + 2*(dalpha[1,0]*fx*fy+dalpha[2,0]*fx*fz+dalpha[2,1]*fy*fz)
        alpha_contrib[i,j] = -5.0e-1*p

        p = 0.0e0
        p = dalphanr[0,0]*fx**2+dalphanr[1,1]*fy**2+dalphanr[2,2]*fz**2
        p = p + 2*(dalphanr[1,0]*fx*fy+dalphanr[2,0]*fx*fz+dalphanr[2,1]*fy*fz)
        alphanr_contrib[i,j] = -5.0e-1*p

        p = 0.0e0
        p = dbeta[0,0,0]*fx**3+dbeta[1,1,1]*fy**3+dbeta[2,2,2]*fz**3
        p = p+3.0e0*(dbeta[0,0,1]*fy*fx**2+dbeta[0,0,2]*fz*fx**2+dbeta[1,0,1]*fx*fy**2)
        p = p+3.0e0*(dbeta[1,1,2]*fz*fy**2+dbeta[2,0,2]*fx*fz**2+dbeta[2,1,2]*fy*fz**2)
        p = p + 6.0e0*dbeta[1,0,2]*fx*fy*fz 
        beta_contrib[i,j] = -(1.0e0/6.0e0)*p

      #Construct each approximation
      if d2_aprox == 1: 
        print(" Aproximation used is only dipole")
        print()
        barrier_2D = nofield_2D+mu_contrib
        barrier_2D = barrier_2D*kcal
      elif d2_aprox == 2:
        print(" Aproximation used is dipole + alpha(el)")
        print()
        barrier_2D = nofield_2D+mu_contrib+alpha_contrib
        barrier_2D = barrier_2D*kcal
      elif d2_aprox == 3:
        print(" Aproximation used is dipole + alpha(el) + beta")
        print()
        barrier_2D = nofield_2D+mu_contrib+alpha_contrib+beta_contrib
        barrier_2D = barrier_2D*kcal
      else : # d2_aprox == 4
        print(" Aproximation used is dipole + alpha(el+nr) + beta")
        print()
        barrier_2D = nofield_2D+mu_contrib+alpha_contrib+alphanr_contrib+beta_contrib
        barrier_2D = barrier_2D*kcal

   elif axis == 2:  
      print(" 2D printing, xz axis")
      print()

      fy = 0.0e0
      for i,fx in enumerate (f_print):
       for j,fz in enumerate (f_print):
        mu_contrib[i,j]      = -1.0e0*(dmu[0]*fx+dmu[1]*fy+dmu[2]*fz)

        p = dalpha[0,0]*fx**2+dalpha[1,1]*fy**2+dalpha[2,2]*fz**2
        p = p + 2*(dalpha[1,0]*fx*fy+dalpha[2,0]*fx*fz+dalpha[2,1]*fy*fz)
        alpha_contrib[i,j] = -5.0e-1*p

        p = 0.0e0
        p = dalphanr[0,0]*fx**2+dalphanr[1,1]*fy**2+dalphanr[2,2]*fz**2
        p = p + 2*(dalphanr[1,0]*fx*fy+dalphanr[2,0]*fx*fz+dalphanr[2,1]*fy*fz)
        alphanr_contrib[i,j] = -5.0e-1*p

        p = 0.0e0
        p = dbeta[0,0,0]*fx**3+dbeta[1,1,1]*fy**3+dbeta[2,2,2]*fz**3
        p = p+3.0e0*(dbeta[0,0,1]*fy*fx**2+dbeta[0,0,2]*fz*fx**2+dbeta[1,0,1]*fx*fy**2)
        p = p+3.0e0*(dbeta[1,1,2]*fz*fy**2+dbeta[2,0,2]*fx*fz**2+dbeta[2,1,2]*fy*fz**2)
        p = p + 6.0e0*dbeta[1,0,2]*fx*fy*fz 
        beta_contrib[i,j] = -(1.0e0/6.0e0)*p

      #Construct each approximation
      if d2_aprox == 1: 
        print(" Aproximation used is only dipole")
        print()
        barrier_2D = nofield_2D+mu_contrib
        barrier_2D = barrier_2D*kcal
      elif d2_aprox == 2:
        print(" Aproximation used is dipole + alpha(el)")
        print()
        barrier_2D = nofield_2D+mu_contrib+alpha_contrib
        barrier_2D = barrier_2D*kcal
      elif d2_aprox == 3:
        print(" Aproximation used is dipole + alpha(el) + beta")
        print()
        barrier_2D = nofield_2D+mu_contrib+alpha_contrib+beta_contrib
        barrier_2D = barrier_2D*kcal
      else : # d2_aprox == 4
        print(" Aproximation used is dipole + alpha(el+nr) + beta")
        print()
        barrier_2D = nofield_2D+mu_contrib+alpha_contrib+alphanr_contrib+beta_contrib
        barrier_2D = barrier_2D*kcal

   elif axis == 3:  
      print(" 2D printing, yz axis")
      print()

      fx = 0.0e0
      for i,fy in enumerate (f_print):
       for j,fz in enumerate (f_print):
        mu_contrib[i,j]      = -1.0e0*(dmu[0]*fx+dmu[1]*fy+dmu[2]*fz)

        p = dalpha[0,0]*fx**2+dalpha[1,1]*fy**2+dalpha[2,2]*fz**2
        p = p + 2*(dalpha[1,0]*fx*fy+dalpha[2,0]*fx*fz+dalpha[2,1]*fy*fz)
        alpha_contrib[i,j] = -5.0e-1*p

        p = 0.0e0
        p = dalphanr[0,0]*fx**2+dalphanr[1,1]*fy**2+dalphanr[2,2]*fz**2
        p = p + 2*(dalphanr[1,0]*fx*fy+dalphanr[2,0]*fx*fz+dalphanr[2,1]*fy*fz)
        alphanr_contrib[i,j] = -5.0e-1*p

        p = 0.0e0
        p = dbeta[0,0,0]*fx**3+dbeta[1,1,1]*fy**3+dbeta[2,2,2]*fz**3
        p = p+3.0e0*(dbeta[0,0,1]*fy*fx**2+dbeta[0,0,2]*fz*fx**2+dbeta[1,0,1]*fx*fy**2)
        p = p+3.0e0*(dbeta[1,1,2]*fz*fy**2+dbeta[2,0,2]*fx*fz**2+dbeta[2,1,2]*fy*fz**2)
        p = p + 6.0e0*dbeta[1,0,2]*fx*fy*fz 
        beta_contrib[i,j] = -(1.0e0/6.0e0)*p

      #Construct each approximation
      if d2_aprox == 1: 
        print(" Aproximation used is only dipole")
        print()
        barrier_2D = nofield_2D+mu_contrib
        barrier_2D = barrier_2D*kcal
      elif d2_aprox == 2:
        print(" Aproximation used is dipole + alpha(el)")
        print()
        barrier_2D = nofield_2D+mu_contrib+alpha_contrib
        barrier_2D = barrier_2D*kcal
      elif d2_aprox == 3:
        print(" Aproximation used is dipole + alpha(el) + beta")
        print()
        barrier_2D = nofield_2D+mu_contrib+alpha_contrib+beta_contrib
        barrier_2D = barrier_2D*kcal
      else : # d2_aprox == 4
        print(" Aproximation used is dipole + alpha(el+nr) + beta")
        print()
        barrier_2D = nofield_2D+mu_contrib+alpha_contrib+alphanr_contrib+beta_contrib
        barrier_2D = barrier_2D*kcal

   else: 
      print(" Bad axis selection in 2D FDB data generation. Stop.")
      print()
      quit()
   
   #Output results for whatever approximation
   #Output small matrix (for scan purposes)
   threshold = 6
   if max(barrier_2D.shape[0],barrier_2D.shape[1]) <= threshold :
     print(np.matrix(barrier_2D))

   else:  #Ouptut big matrix (element by element printing)
     #for each first label, print all 2nd label.
     #eg. for xy calculation prints: for x1 all y, then for x2 all y... up to for x_last prints all y
     for i in range(fpoints):
      for j in range(fpoints):
        print(barrier_2D[i,j])

else:
   print()
   print(" Bad input for FDB approximation in 2D data generation. Stop. ")
   print()
   file1.close()
   quit()

file1.close()

### End of 2D FDB data generation ###
### End of script ###


