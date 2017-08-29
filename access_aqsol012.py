# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 19:18:12 2016

@author: PC
"""
import sys
import os
if os.environ['COMPUTERNAME'] == 'WAHLINJOHANSEN':
    sys.path.append ("E:\\NTNU\\python_lib")
else:
    sys.path.append ("C:\Work\python_lib\solutionCalcs")
import ctypes
import numpy as np
#import iapws

if os.environ['COMPUTERNAME'] == 'WAHLINJOHANSEN':
    print(os.environ['COMPUTERNAME'])
    mydll=ctypes.WinDLL("E:\\NTNU\\python_lib\\AQSOL012-64.dll")
elif os.environ['COMPUTERNAME'] == 'IBATMBWAHLIN':
    mydll=ctypes.WinDLL("C:\\Users\\wahlin\\Documents\\Vegvesen\\python_lib\\solutionCalcs\\AQSOL012.dll")
else:
    mydll=ctypes.WinDLL("C:\\python_scripts\\solutionCalcs\\AQSOL012-64.dll")

M_Na_plus= 22.98977
M_K_plus = 39.09831
M_Mg_plus = 24.30506
M_Ca_plus = 40.0784
M_Cl_minus = 35.4532

def get_aw(mol_nacl, mol_mgcl, mol_cacl, mol_kcl,mass_water, temp):
    gramH2O = mass_water
    gramNa = mol_nacl*M_Na_plus
    gramK = mol_kcl*M_K_plus
    gramMg = mol_mgcl*M_Mg_plus
    gramCa = mol_cacl*M_Ca_plus
    gramH = 0
    gramCl = mol_nacl*M_Cl_minus + 2*mol_mgcl*M_Cl_minus + 2*mol_cacl*M_Cl_minus + mol_kcl*M_Cl_minus
    gramSO4 = 0
    gramHSO4 = 0
    gramOH = 0

    #mydll=ctypes.WinDLL("C:\\Users\\wahlin\\Documents\\Work\\python_lib\\AQSOL012.dll")
    
    nk = 12
    nsalts = 69
    nsalte = nsalts - nk
    nsalt2 = 2 * nsalte
    

    SI = np.zeros(nsalt2, dtype=ctypes.c_double) 
    gram = np.zeros(nk,dtype=ctypes.c_double) 
    eqgram = np.zeros(nk,dtype=ctypes.c_double) 
    fgram = np.zeros(nsalts,dtype=ctypes.c_double) 
    Faktor = np.zeros(nsalts,dtype=ctypes.c_double) 
    ENTL = np.zeros(1,dtype=ctypes.c_double) 
    ENTS = np.zeros(1, dtype=ctypes.c_double) 
    PBUB = np.zeros(1, dtype=ctypes.c_double) 
    ie = np.zeros(1,dtype=ctypes.c_long)
    CPL = np.zeros(1, dtype=ctypes.c_double) 
    pH = np.zeros(1, dtype=ctypes.c_double) 
    pHC = np.zeros(1, dtype=ctypes.c_double) 
    AW = np.zeros(1, dtype=ctypes.c_double) 
    IStr = np.zeros(1, dtype=ctypes.c_double) 
    ie = np.zeros(1, dtype=ctypes.c_long) 
    temperatur = np.zeros(1, dtype=ctypes.c_double) 
    timeout = np.zeros(1, dtype=ctypes.c_double) 
    
    temperatur[0]=temp
    Faktor[0:57]=1
    timeout[0] = 25
    
    
    gram[0]=gramH2O
    gram[1]=gramNa
    gram[2]=gramK
    gram[3]=gramMg
    gram[4]=gramCa
    gram[5]=gramH
    gram[6]=0
    gram[7]=gramCl
    gram[8]=gramSO4
    gram[9]=gramHSO4
    gram[10]=gramOH
    gram[11]=0
    
    
    mydll.EQLB(temperatur.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),gram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), Faktor.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),eqgram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), fgram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), SI.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ENTL.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ENTS.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), CPL.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),PBUB.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),pH.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),pHC.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),AW.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),IStr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),timeout.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),ie.ctypes.data_as(ctypes.POINTER(ctypes.c_long)))
    
    if ie ==0:
        activity=AW
    else:
        activity = np.nan
    
    return (activity)
    
def get_solid_ice(mol_nacl, mol_mgcl, mol_cacl, mol_kcl,mass_water, temp):
    gramH2O = mass_water
    gramNa = mol_nacl*M_Na_plus
    gramK = mol_kcl*M_K_plus
    gramMg = mol_mgcl*M_Mg_plus
    gramCa = mol_cacl*M_Ca_plus
    gramH = 0
    gramCl = mol_nacl*M_Cl_minus + 2*mol_mgcl*M_Cl_minus + 2*mol_cacl*M_Cl_minus + mol_kcl*M_Cl_minus
    gramSO4 = 0
    gramHSO4 = 0
    gramOH = 0

    #mydll=ctypes.WinDLL("C:\\Users\\wahlin\\Documents\\Work\\python_lib\\AQSOL012.dll")
    
    nk = 12
    nsalts = 69
    nsalte = nsalts - nk
    nsalt2 = 2 * nsalte
    

    SI = np.zeros(nsalt2, dtype=ctypes.c_double) 
    gram = np.zeros(nk,dtype=ctypes.c_double) 
    eqgram = np.zeros(nk,dtype=ctypes.c_double) 
    fgram = np.zeros(nsalts,dtype=ctypes.c_double) 
    Faktor = np.zeros(nsalts,dtype=ctypes.c_double) 
    ENTL = np.zeros(1,dtype=ctypes.c_double) 
    ENTS = np.zeros(1, dtype=ctypes.c_double) 
    PBUB = np.zeros(1, dtype=ctypes.c_double) 
    ie = np.zeros(1,dtype=ctypes.c_long)
    CPL = np.zeros(1, dtype=ctypes.c_double) 
    pH = np.zeros(1, dtype=ctypes.c_double) 
    pHC = np.zeros(1, dtype=ctypes.c_double) 
    AW = np.zeros(1, dtype=ctypes.c_double) 
    IStr = np.zeros(1, dtype=ctypes.c_double) 
    ie = np.zeros(1, dtype=ctypes.c_long) 
    temperatur = np.zeros(1, dtype=ctypes.c_double) 
    timeout = np.zeros(1, dtype=ctypes.c_double) 
    
    temperatur[0]=temp
    Faktor[0:57]=1
    timeout[0] = 25
    
    
    gram[0]=gramH2O
    gram[1]=gramNa
    gram[2]=gramK
    gram[3]=gramMg
    gram[4]=gramCa
    gram[5]=gramH
    gram[6]=0
    gram[7]=gramCl
    gram[8]=gramSO4
    gram[9]=gramHSO4
    gram[10]=gramOH
    gram[11]=0
    
   
    mydll.EQLB(temperatur.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),gram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), Faktor.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),eqgram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), fgram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), SI.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ENTL.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ENTS.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), CPL.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),PBUB.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),pH.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),pHC.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),AW.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),IStr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),timeout.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),ie.ctypes.data_as(ctypes.POINTER(ctypes.c_long)))
    mass_ice = fgram[-4]
    return (mass_ice)



    
def get_SI_ice(mol_nacl, mol_mgcl, mol_cacl, mol_kcl,mass_water, temp):
    
#    mol_nacl = 2
#    mol_mgcl = 0
#    mol_cacl= 0
#    mol_kcl = 0
#    mass_water =1000
#    temp = -14
##    
    gramH2O = mass_water
    gramNa = mol_nacl*M_Na_plus
    gramK = mol_kcl*M_K_plus
    gramMg = mol_mgcl*M_Mg_plus
    gramCa = mol_cacl*M_Ca_plus
    gramH = 0
    gramCl = mol_nacl*M_Cl_minus + 2*mol_mgcl*M_Cl_minus + 2*mol_cacl*M_Cl_minus + mol_kcl*M_Cl_minus
    gramSO4 = 0
    gramHSO4 = 0
    gramOH = 0
    
    
    #mydll=ctypes.WinDLL("C:\\Users\\wahlin\\Documents\\Work\\python_lib\\AQSOL012.dll")
    
    nk = 12
    nsalts = 69
    #resl = 3 * nsalts + 14 - nk
    nsalte = nsalts - nk
    nsalt2 = 2 * nsalte
    
    
   # result = np.zeros(resl, dtype=ctypes.c_double) 
    SI = np.zeros(nsalt2, dtype=ctypes.c_double) 
    gram = np.zeros(nk,dtype=ctypes.c_double) 
    eqgram = np.zeros(nk,dtype=ctypes.c_double) 
    fgram = np.zeros(nsalts,dtype=ctypes.c_double) 
    Faktor = np.zeros(nsalts,dtype=ctypes.c_double) 
    ENTL = np.zeros(1,dtype=ctypes.c_double) 
    ENTS = np.zeros(1, dtype=ctypes.c_double) 
    PBUB = np.zeros(1, dtype=ctypes.c_double) 
    ie = np.zeros(1,dtype=ctypes.c_long)
    CPL = np.zeros(1, dtype=ctypes.c_double) 
    pH = np.zeros(1, dtype=ctypes.c_double) 
    pHC = np.zeros(1, dtype=ctypes.c_double) 
    AW = np.zeros(1, dtype=ctypes.c_double) 
    IStr = np.zeros(1, dtype=ctypes.c_double) 
    ie = np.zeros(1, dtype=ctypes.c_long) 
    temperatur = np.zeros(1, dtype=ctypes.c_double) 
    timeout = np.zeros(1, dtype=ctypes.c_double) 
    
    temperatur[0]=temp
    Faktor[0:57]=1
    timeout[0] = 25
    
    
    gram[0]=gramH2O
    gram[1]=gramNa
    gram[2]=gramK
    gram[3]=gramMg
    gram[4]=gramCa
    gram[5]=gramH
    gram[6]=0
    gram[7]=gramCl
    gram[8]=gramSO4
    gram[9]=gramHSO4
    gram[10]=gramOH
    gram[11]=0
    
    
    mydll.EQLB(temperatur.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),gram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), Faktor.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),eqgram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), fgram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), SI.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ENTL.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ENTS.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), CPL.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),PBUB.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),pH.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),pHC.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),AW.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),IStr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),timeout.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),ie.ctypes.data_as(ctypes.POINTER(ctypes.c_long)))
    
    
    if ie ==0:
        SI_ice = SI[-4]
    else:
        SI_ice = np.nan
    
    return (SI_ice)
    
    
def get_alldata(mol_nacl, mol_mgcl, mol_cacl, mol_kcl,mass_water, temp):

   
    gramH2O = mass_water
    gramNa = mol_nacl*M_Na_plus
    gramK = mol_kcl*M_K_plus
    gramMg = mol_mgcl*M_Mg_plus
    gramCa = mol_cacl*M_Ca_plus
    gramH = 0
    gramCl = mol_nacl*M_Cl_minus + 2*mol_mgcl*M_Cl_minus + 2*mol_cacl*M_Cl_minus + mol_kcl*M_Cl_minus
    gramSO4 = 0
    gramHSO4 = 0
    gramOH = 0
    
    
    #mydll=ctypes.WinDLL("C:\\Users\\wahlin\\Documents\\Work\\python_lib\\AQSOL012.dll")
    
    nk = 12
    nsalts = 69
    #resl = 3 * nsalts + 14 - nk
    nsalte = nsalts - nk
    nsalt2 = 2 * nsalte
    
    
    
   # result = np.zeros(resl, dtype=ctypes.c_double) 
    SI = np.zeros(nsalt2, dtype=ctypes.c_double) 
    gram = np.zeros(nk,dtype=ctypes.c_double) 
    eqgram = np.zeros(nk,dtype=ctypes.c_double) 
    fgram = np.zeros(nsalts,dtype=ctypes.c_double) 
    Faktor = np.zeros(nsalts,dtype=ctypes.c_double) 
    ENTL = np.zeros(1,dtype=ctypes.c_double) 
    ENTS = np.zeros(1, dtype=ctypes.c_double) 
    PBUB = np.zeros(1, dtype=ctypes.c_double) 
    ie = np.zeros(1,dtype=ctypes.c_long)
    CPL = np.zeros(1, dtype=ctypes.c_double) 
    pH = np.zeros(1, dtype=ctypes.c_double) 
    pHC = np.zeros(1, dtype=ctypes.c_double) 
    AW = np.zeros(1, dtype=ctypes.c_double) 
    IStr = np.zeros(1, dtype=ctypes.c_double) 
    ie = np.zeros(1, dtype=ctypes.c_long) 
    temperatur = np.zeros(1, dtype=ctypes.c_double) 
    timeout = np.zeros(1, dtype=ctypes.c_double) 
    
    temperatur[0]=temp
    Faktor[0:57]=1
    timeout[0] = 25
    
    
    gram[0]=gramH2O
    gram[1]=gramNa
    gram[2]=gramK
    gram[3]=gramMg
    gram[4]=gramCa
    gram[5]=gramH
    gram[6]=0
    gram[7]=gramCl
    gram[8]=gramSO4
    gram[9]=gramHSO4
    gram[10]=gramOH
    gram[11]=0
    
    
    
    mydll.EQLB(temperatur.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),gram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), Faktor.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),eqgram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), fgram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), SI.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ENTL.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ENTS.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), CPL.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),PBUB.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),pH.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),pHC.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),AW.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),IStr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),timeout.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),ie.ctypes.data_as(ctypes.POINTER(ctypes.c_long)))
    
    titles_properties = ['Dissolved species','Solid phases','Saturation index before equlibration',
                     'Saturation index after equlibration','Enhtalpy of liquid', 'Enthalpy of solid',
                     'Heat capacity liquid','Bubble pressure','pH','Mean ionic molal activity coefficient',
                     'Ionic strength','Water activity','Error message']
    labels_dissolved = ['Water', "Na+", "K+", "Mg++","Ca++", "H+","Cl-", "SO4--","HSO4-","OH-"]
    labels_solids = ["NaCl","NaCl*2H2O","Na2SO4","Na2SO4*10H2O","NaHSO4","NaOH","NaOH*H2O","2NaOH*7H2O","KCl ","K2SO4",
                      "KHSO4","K2SO4*KHSO4","KOH","KOH*H2O","KOH*2H2O","NaK3(SO4)2","MgCl2","MgCl2*2H2O ","MgCl2*4H2O",
                      "MgCl2*6H2O","MgCl2*8H2O","MgCl2*12H2O","KCl*MgCl2*6H2O","MgSO4","MgSO4*H2O","MgSO4*6H2O","MgSO4*7H2O",
                      "MgSO4*12H2O","Na2SO4*MgSO4*4H2O","3Na2SO4*MgSO4","2Na2SO4*2MgSO4*5H2O","KCl*MgSO4*3H2O","K2SO4*2MgSO4",
                      "K2SO4*MgSO4*4H2O","K2SO4*MgSO4*6H2O","Mg(OH)2","CaCl2","CaCl2*H2O","CaCl2*2H2O","CaCl2*4H2O","CaCl2*6H2O",
                      "KCl*CaCl2","2MgCl2*CaCl2*12H2O","MgCl2*2CaCl2*6H2O","CaSO4","2CaSO4*H2O","CaSO4*2H2O","Na2SO4*CaSO4",
                      "K2SO4*MgSO4*2CaSO4*2H2O","K2SO4*CaSO4*H2O","K2SO4*5CaSO4*H2O","Ca(OH)2","CaCl2*Ca(OH)2*H2O","Ice",
                      "H2SO4","HCl(g)","H2O(g)"]
    
    data_dict={}
    
    data_dict[titles_properties[0]]= dict(zip(labels_dissolved,eqgram[[0,1,2,3,4,5,7,8,9,10]]))
    data_dict[titles_properties[1]]= dict(zip(labels_solids,fgram[12:]))
    data_dict[titles_properties[2]]= dict(zip(labels_solids,SI[1:58]))
    data_dict[titles_properties[3]]=dict(zip(labels_solids,SI[58:]))
    data_dict[titles_properties[4]]= ENTL[0]
    data_dict[titles_properties[5]]= ENTS[0]
    data_dict[titles_properties[6]]= CPL[0]
    data_dict[titles_properties[7]]= PBUB[0]
    data_dict[titles_properties[8]]= pH[0]
    data_dict[titles_properties[9]]= pHC[0]
    data_dict[titles_properties[10]]= IStr[0]
    data_dict[titles_properties[11]]= AW[0]
    data_dict[titles_properties[12]]= ie[0]
    
    
    return data_dict

def get_sold_mass(mol_nacl, mol_mgcl, mol_cacl, mol_kcl,mass_water, temp):

   
    gramH2O = mass_water
    gramNa = mol_nacl*M_Na_plus
    gramK = mol_kcl*M_K_plus
    gramMg = mol_mgcl*M_Mg_plus
    gramCa = mol_cacl*M_Ca_plus
    gramH = 0
    gramCl = mol_nacl*M_Cl_minus + 2*mol_mgcl*M_Cl_minus + 2*mol_cacl*M_Cl_minus + mol_kcl*M_Cl_minus
    gramSO4 = 0
    gramHSO4 = 0
    gramOH = 0
    
    
    #mydll=ctypes.WinDLL("C:\\Users\\wahlin\\Documents\\Work\\python_lib\\AQSOL012.dll")
    
    nk = 12
    nsalts = 69
    #resl = 3 * nsalts + 14 - nk
    nsalte = nsalts - nk
    nsalt2 = 2 * nsalte
    
    
    
   # result = np.zeros(resl, dtype=ctypes.c_double) 
    SI = np.zeros(nsalt2, dtype=ctypes.c_double) 
    gram = np.zeros(nk,dtype=ctypes.c_double) 
    eqgram = np.zeros(nk,dtype=ctypes.c_double) 
    fgram = np.zeros(nsalts,dtype=ctypes.c_double) 
    Faktor = np.zeros(nsalts,dtype=ctypes.c_double) 
    ENTL = np.zeros(1,dtype=ctypes.c_double) 
    ENTS = np.zeros(1, dtype=ctypes.c_double) 
    PBUB = np.zeros(1, dtype=ctypes.c_double) 
    ie = np.zeros(1,dtype=ctypes.c_long)
    CPL = np.zeros(1, dtype=ctypes.c_double) 
    pH = np.zeros(1, dtype=ctypes.c_double) 
    pHC = np.zeros(1, dtype=ctypes.c_double) 
    AW = np.zeros(1, dtype=ctypes.c_double) 
    IStr = np.zeros(1, dtype=ctypes.c_double) 
    ie = np.zeros(1, dtype=ctypes.c_long) 
    temperatur = np.zeros(1, dtype=ctypes.c_double) 
    timeout = np.zeros(1, dtype=ctypes.c_double) 
    
    temperatur[0]=temp
    Faktor[0:57]=1
    timeout[0] = 25
    
    
    gram[0]=gramH2O
    gram[1]=gramNa
    gram[2]=gramK
    gram[3]=gramMg
    gram[4]=gramCa
    gram[5]=gramH
    gram[6]=0
    gram[7]=gramCl
    gram[8]=gramSO4
    gram[9]=gramHSO4
    gram[10]=gramOH
    gram[11]=0
    
    
    
    mydll.EQLB(temperatur.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),gram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), Faktor.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),eqgram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), fgram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), SI.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ENTL.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ENTS.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), CPL.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),PBUB.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),pH.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),pHC.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),AW.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),IStr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),timeout.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),ie.ctypes.data_as(ctypes.POINTER(ctypes.c_long)))
    
    
    
    return (fgram)
    
    
def get_thdyn_properties(mol_nacl, mol_mgcl, mol_cacl, mol_kcl,mass_water, temp):
    gramH2O = mass_water
    gramNa = mol_nacl*M_Na_plus
    gramK = mol_kcl*M_K_plus
    gramMg = mol_mgcl*M_Mg_plus
    gramCa = mol_cacl*M_Ca_plus
    gramH = 0
    gramCl = mol_nacl*M_Cl_minus + 2*mol_mgcl*M_Cl_minus + 2*mol_cacl*M_Cl_minus + mol_kcl*M_Cl_minus
    gramSO4 = 0
    gramHSO4 = 0
    gramOH = 0

    #mydll=ctypes.WinDLL("C:\\Users\\wahlin\\Documents\\Work\\python_lib\\AQSOL012.dll")
    
    nk = 12
    nsalts = 69
    nsalte = nsalts - nk
    nsalt2 = 2 * nsalte
    

    SI = np.zeros(nsalt2, dtype=ctypes.c_double) 
    gram = np.zeros(nk,dtype=ctypes.c_double) 
    eqgram = np.zeros(nk,dtype=ctypes.c_double) 
    fgram = np.zeros(nsalts,dtype=ctypes.c_double) 
    Faktor = np.zeros(nsalts,dtype=ctypes.c_double) 
    ENTL = np.zeros(1,dtype=ctypes.c_double) 
    ENTS = np.zeros(1, dtype=ctypes.c_double) 
    PBUB = np.zeros(1, dtype=ctypes.c_double) 
    ie = np.zeros(1,dtype=ctypes.c_long)
    CPL = np.zeros(1, dtype=ctypes.c_double) 
    pH = np.zeros(1, dtype=ctypes.c_double) 
    pHC = np.zeros(1, dtype=ctypes.c_double) 
    AW = np.zeros(1, dtype=ctypes.c_double) 
    IStr = np.zeros(1, dtype=ctypes.c_double) 
    ie = np.zeros(1, dtype=ctypes.c_long) 
    temperatur = np.zeros(1, dtype=ctypes.c_double) 
    timeout = np.zeros(1, dtype=ctypes.c_double) 
    
    temperatur[0]=temp
    Faktor[0:57]=1
    timeout[0] = 25
    
    
    gram[0]=gramH2O
    gram[1]=gramNa
    gram[2]=gramK
    gram[3]=gramMg
    gram[4]=gramCa
    gram[5]=gramH
    gram[6]=0
    gram[7]=gramCl
    gram[8]=gramSO4
    gram[9]=gramHSO4
    gram[10]=gramOH
    gram[11]=0
    
    
    mydll.EQLB(temperatur.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),gram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), Faktor.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),eqgram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), fgram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), SI.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ENTL.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ENTS.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), CPL.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),PBUB.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),pH.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),pHC.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),AW.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),IStr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),timeout.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),ie.ctypes.data_as(ctypes.POINTER(ctypes.c_long)))
    
    return (ENTL, ENTS, CPL)