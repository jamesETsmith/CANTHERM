#!/usr/bin/env python

import sys
sys.path.append('/home/sandeeps/site-packages')
sys.path.append('/home/sandeeps/site-packages/Numeric')
import readGeomFc
import pdb
import math
from numpy import *
from scipy import *

class CanTherm:
 CalcType = ''
 ReacType = ''
 Temp = []
 MoleculeList = []
 Entropy = []
 Thermal = []
 Cp = []
 Parition = []
 scale = 0.0
 #CBSQB3 E for H, N, O, C, P
 atomEcbsqb3 = {'H':-0.499818 , 'N':-54.520543 , 'O':-74.987624 , 'C':-37.785385 , 'P':-340.817186}
 #G3 E for H, N, O, C, P
 atomEg3 = {'H':-0.5010030, 'N':-54.564343, 'O':-75.030991, 'C':-37.827717, 'P':-341.116432}
 #expt H contains H + TC + SOC (spin orbital correction)
 atomH = {'H':50.62 , 'N':111.49 , 'O':58.163 , 'C':169.8147 }
 #BAC for C-H C-C C=C C.TB.C  O-H  C-O C=O
 bondC = [-0.11, -0.3, -0.08, -0.64, 0.02, 0.33, 0.55]



def main():
  data = CanTherm()
  inputFile = open(sys.argv[1],'r') 
  oFile = open('cantherm.out','w')
  readGeomFc.readInputFile(inputFile,data)

  data.Entropy=len(data.MoleculeList)*len(data.Temp)*[0.0]
  data.Cp=len(data.MoleculeList)*len(data.Temp)*[0.0]
  data.Thermal=len(data.MoleculeList)*len(data.Temp)*[0.0]
  data.Partition=len(data.MoleculeList)*len(data.Temp)*[1.0]
  Entropy = data.Entropy
  Cp = data.Cp
  Thermal = data.Thermal
  Partition = data.Partition

  for i in range(len(data.MoleculeList)):
     molecule = data.MoleculeList[i]
     oFile.write('Molecule '+str(i+1)+':\n')
     oFile.write('-----------\n\n')
     molecule.printData(oFile)

     oFile.write('\nThermodynamic Data\n')

     Temp = data.Temp
     #translation
     (ent,cp,dh) = molecule.getTranslationThermo(oFile,data.Temp) 
     for j in range(len(Temp)):
         Entropy[i*len(Temp)+j]=Entropy[i*len(Temp)+j]+ent[j]
         Cp[i*len(Temp)+j]=Cp[i*len(Temp)+j]+cp[j]
#         Thermal[i*len(Temp)+j]=Thermal[i*len(Temp)+j]+dh[j]
  
     #vibrational
     (ent,cp,dh,q) = molecule.getVibrationalThermo(oFile,data.Temp,data.scale) 
     for j in range(len(Temp)):
         Entropy[i*len(Temp)+j]=Entropy[i*len(Temp)+j]+ent[j]
         Cp[i*len(Temp)+j]=Cp[i*len(Temp)+j]+cp[j]
         Thermal[i*len(Temp)+j]=Thermal[i*len(Temp)+j]+dh[j]
         #Partition[i*len(Temp)+j] = Partition[i*len(Temp)+j]*q[j]
         #print '%12.2f'%float(ent[j]),
     #print '\n'

     #Internal rotational
     if molecule.numRotors != 0:
      (ent,cp,dh,q) = molecule.getIntRotationalThermo_Q(oFile,data.Temp) 
      for j in range(len(Temp)):
         Entropy[i*len(Temp)+j]=Entropy[i*len(Temp)+j]+ent[j]
         Cp[i*len(Temp)+j]=Cp[i*len(Temp)+j]+cp[j]
         Thermal[i*len(Temp)+j]=Thermal[i*len(Temp)+j]+dh[j]
         #Partition[i*len(Temp)+j] = Partition[i*len(Temp)+j]*q[j]
         #print '%12.2f'%float(ent[j]),
     #print '\n'

     #External rotational
     (ent,cp,dh) = molecule.getExtRotationalThermo(oFile,data.Temp) 
     for j in range(len(Temp)):
         Entropy[i*len(Temp)+j]=Entropy[i*len(Temp)+j]+ent[j]
         Cp[i*len(Temp)+j]=Cp[i*len(Temp)+j]+cp[j]
         Thermal[i*len(Temp)+j]=Thermal[i*len(Temp)+j]+dh[j]
#         Partition[i*len(Temp)+j] = Thermal[i*len(Temp)+j]+q[j]

     for j in range(len(Temp)):
         Entropy[i*len(Temp)+j]=Entropy[i*len(Temp)+j]+1.985*math.log(molecule.nelec)

     #print Enthalpy
     H = molecule.Energy
     atoms = readGeomFc.getAtoms(molecule.Mass)
     atomsH = 0.0
     if molecule.Etype == 'cbsqb3':
        atomE = data.atomEcbsqb3
     if molecule.Etype == 'g3':
        atomE = data.atomEg3
     for atom in atoms:
         H -= atomE[atom]
         atomsH += data.atomH[atom]
     H = H*627.5095+atomsH

     if molecule.Etype == 'cbsqb3':
       b = 0
       for bonds in molecule.bonds:
         H += bonds*data.bondC[b]
         b += 1

     H += Thermal[i*len(Temp)+0]
     print '%12.2f'%H + '%12.2f'%Entropy[i*len(Temp)+0],
     for c in range(1,8):
        print '%12.2f'%Cp[i*len(Temp)+c],
     print '\n'

     #for c in range(len(Temp)):
     #   print '%12.2e'%Partition[i*len(Temp)+c],
     #print

  if len(data.MoleculeList) == 1:
     return

  #fit the rate coefficient
  A = matrix(zeros((len(Temp),3),dtype=float))
  y = matrix(zeros((len(Temp),1),dtype=float))

  rate = [0.0]*len(Temp)
  for j in range(len(Temp)):
    if (data.ReacType == 'Unimol'):
      rate[j] = (1.381e-23*Temp[j]/6.626e-34)*math.exp((Entropy[len(Temp)+j]-Entropy[j])/1.985)*math.exp(-(data.MoleculeList[1].Energy - data.MoleculeList[0].Energy)*627.5095*1.0e3/1.985/Temp[j])
      #wigner correction
      rate[j] *= (1.0 +1.0/24.0 * (1.44*data.MoleculeList[1].imagFreq/Temp[j]) )**2
      A[j,:] = mat([1.0, math.log(Temp[j]), -1.0/1.985/Temp[j]])
      y[j] = log(rate[j])
  b = linalg.inv(transpose(A)*A)*(transpose(A)*y)
  oFile.write('\n\nRate Data\n')
  oFile.write('r = A*(T/1000)^n*exp(-Ea/R/T)'+'%12.2e'%(exp(b[0])*1000.0**float(b[1]))+'%10.2f'%b[1]+'%12.2f'%(b[2]/1.0e3)+
'\n')
  oFile.write('r = A*T^n*exp(-Ea/R/T)'+'%12.2e'%(exp(b[0]))+'%10.2f'%b[1]+'%12.2f'%(b[2]/1.0e3)+'\n')
  oFile.write('%12s'%'Temperature'+'%12s\n'%'Rate')
  for j in range(len(Temp)):
      fitrate = exp(b[0])*Temp[j]**float(b[1])*exp(-b[2]/1.985/Temp[j])
      oFile.write('%12.2f'%Temp[j]+'%12.2e'%rate[j]+'%12.2e\n'%fitrate)
  oFile.write('\n\n')
  oFile.close()
  

if __name__ == "__main__":
   main()

