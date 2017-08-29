# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 10:29:07 2016

@author: wahlin
"""
import sys
sys.path.append ("C:/Users/wahlin/Documents/Work/python_lib")
sys.path.append ("E:/NTNU/python_lib")

import scipy.constants
import numpy as np
from iapws import IAPWS95
from iapws import _Ice
import access_aqsol012 as aq12
from imp import reload
#import matplotlib.pyplot as plt

molar_w_water = 18.01533 # g/mol
mol_kg_water = 1000/molar_w_water


def estimate_freezpoint_aw(aw):
    water = IAPWS95(T=273.15, P=0.101325)
    ice=_Ice(T=273.15,P=0.101325)

    enth_fus = 1000*(water.h-ice["h"])/mol_kg_water # enthalphy of fusion for water at 0 degrees centigrade
    cP_diff = 1000*(water.cp-ice["cp"])/mol_kg_water # heat capacity difference between ice and water at 0 centigrade
    R = scipy.constants.R
    T_f = 273.15
    
    fpd = -(enth_fus - 2*R*T_f*np.log(aw)- np.sqrt(2*cP_diff*T_f**2*R*np.log(aw)+enth_fus**2))/(2*(enth_fus/T_f+cP_diff/2-R*np.log(aw)))
    
    return fpd
    
def get_vapour_pressure(temp, state):
    # Accepts a temperature in degrees Celcius and and state ("ice" or "water") and returns the activity of that state at that temperature (Following Murphy and Koop (2005))    
    TK = temp+273.15
    if state=="water":
        p_vap = np.exp(54.842763 - 6763.22/TK -4.210*np.log(TK) + 0.000367*TK + np.tanh(0.0415*(TK-218.8))*(53.878-1331.22/TK-9.44523*np.log(TK)+0.014025*TK)) # From Murphy and Koop (2005). Valid range: 123 < T < 332 K. Source is reliable and can be trusted
    elif state=="ice":
        p_vap = np.exp(9.550426 - 5723.256/TK + 3.53068*np.log(TK)-0.00728332*TK) # From Murphy and Koop (2005). Valid range: >110 K. Source is reliable and can be trusted.
    else:
        print("Unknown state, please choose between water and ice.")
        p_vap=np.nan
    return p_vap
    
def estimate_aw_from_fp(fp,temp):
    # Accepts temperatures in degrees C (freezing porint and evaluation temperature) and returns the approximate solution activity
        
    dT = np.abs(fp)
    temp = 273.15+temp
    Tf=273.15
    water = IAPWS95(T=temp, P=0.101325)
    ice=_Ice(T=temp,P=0.101325)

    enth_fus = 1000*(water.h-ice["h"]) # enthalphy of fusion for water at 0 degrees centigrade
    cP_diff = 1000*(water.cp-ice["cp"]) # heat capacity difference between ice and water at 0 centigrade
    R = scipy.constants.R
  
  
    term1=enth_fus*(1/Tf-1/(Tf-dT))
    term2=cP_diff*(np.log((Tf-dT)/Tf)-dT/(Tf-dT))
    t12 = (term1+term2)/(R*temp)
    aw=np.exp(t12)
    return aw  
    
def calc_icemelt_cap(n_nacl, n_mgcl, n_cacl, n_kcl, m_water, temp):
    ''' Composition is a 4x1 vector giving the molar concentration of
NaCl, MgCl, CaCl, KCl in the given order. Temp is simply the temperature. The script uses 1000 gram of water 
into which the salts are mixed. This should be adjusted when I come up with a good idea of what is most relevant '''
    composition=np.zeros(4)
    composition[0]= n_nacl
    composition[1] = n_mgcl
    composition[2]= n_cacl
    composition[3] = n_kcl
    
    m_melt = float(0)
    SI_ice = aq12.get_SI_ice(composition[0],composition[1],composition[2],composition[3],m_water+m_melt,temp)
    if np.isnan(SI_ice):
        while np.isnan(SI_ice):
            print(SI_ice)
            m_melt=m_melt+0.1
            SI_ice = aq12.get_SI_ice(composition[0],composition[1],composition[2],composition[3],m_water+m_melt,temp)
            print(SI_ice)

    if np.abs(SI_ice-1) > 0.001:
            while np.abs(SI_ice-1) > 0.001:
                m_melt=m_melt+10000
                SI_ice = aq12.get_SI_ice(composition[0],composition[1],composition[2],composition[3],m_water+m_melt,temp)
            m_melt = m_melt -10000   
            SI_ice = aq12.get_SI_ice(composition[0],composition[1],composition[2],composition[3],m_water+m_melt,temp)        
        
            while np.abs(SI_ice-1) > 0.001:
                m_melt=m_melt+1000
                SI_ice = aq12.get_SI_ice(composition[0],composition[1],composition[2],composition[3],m_water+m_melt,temp)
            m_melt = m_melt -1000   
            SI_ice = aq12.get_SI_ice(composition[0],composition[1],composition[2],composition[3],m_water+m_melt,temp)
            
            while np.abs(SI_ice-1) > 0.001:
                m_melt=m_melt+100
                SI_ice = aq12.get_SI_ice(composition[0],composition[1],composition[2],composition[3],m_water+m_melt,temp)
            m_melt = m_melt -100   
            SI_ice = aq12.get_SI_ice(composition[0],composition[1],composition[2],composition[3],m_water+m_melt,temp)

            while np.abs(SI_ice-1) > 0.001:
                m_melt=m_melt+10
                SI_ice = aq12.get_SI_ice(composition[0],composition[1],composition[2],composition[3],m_water+m_melt,temp)
                
            m_melt = m_melt -10  
            SI_ice = aq12.get_SI_ice(composition[0],composition[1],composition[2],composition[3],m_water+m_melt,temp)

            while np.abs(SI_ice-1) > 0.001:
                m_melt=m_melt+1
                SI_ice = aq12.get_SI_ice(composition[0],composition[1],composition[2],composition[3],m_water+m_melt,temp)
                
            m_melt = m_melt -1 
            SI_ice = aq12.get_SI_ice(composition[0],composition[1],composition[2],composition[3],m_water+m_melt,temp)

            while np.abs(SI_ice-1) > 0.001:
                m_melt=m_melt+0.1
                SI_ice = aq12.get_SI_ice(composition[0],composition[1],composition[2],composition[3],m_water+m_melt,temp)
                
            m_melt = m_melt -0.1 
            SI_ice = aq12.get_SI_ice(composition[0],composition[1],composition[2],composition[3],m_water+m_melt,temp)
            
            while np.abs(SI_ice-1) > 0.001:
                m_melt=m_melt+0.01
                SI_ice = aq12.get_SI_ice(composition[0],composition[1],composition[2],composition[3],m_water+m_melt,temp)
                
            m_melt = m_melt -0.01 
            SI_ice = aq12.get_SI_ice(composition[0],composition[1],composition[2],composition[3],m_water+m_melt,temp)
        
            melt_cap=m_melt
    else:
        melt_cap=0
        
    return (melt_cap)
    
    
def find_solubility(n_nacl, n_mgcl, n_cacl, n_kcl,temp):
    ''' Composition is 4x1 vectors giving the molar concentration of
NaCl, MgCl, CaCl, KCl in the given order. Temp is simply the temperature. The script begins with the given composition dissolved in a 1000 gram of water. 
If the given composition does not fully dissolve in 1000g of water, the function returns nan. This should be adjusted when I come up with a good idea of what is most relevant '''
    
    composition = [n_nacl, n_mgcl, n_cacl, n_kcl]
    m_evap = 0
    prod =aq12.get_alldata(composition[0],composition[1],composition[2],composition[3],1000,temp)
    print(prod['Error message'])
    if prod['Error message'] == 0:    
        prod['Solid phases']['Ice']=0
        tot_solid=sum(prod['Solid phases'].values())
    else:
        tot_solid =10
    #print(tot_solid)
    if np.abs(tot_solid) < 0.001:       
            while np.abs(tot_solid) < 0.001:
                m_evap=m_evap+100
                prod =aq12.get_alldata(composition[0],composition[1],composition[2],composition[3],1000-m_evap,temp)
                if prod['Error message'] == 0:    
                    prod['Solid phases']['Ice']=0
                    tot_solid=sum(prod['Solid phases'].values())
                else:
                    tot_solid =10
            m_evap = m_evap -100   
            prod =aq12.get_alldata(composition[0],composition[1],composition[2],composition[3],1000-m_evap,temp)
            prod['Solid phases']['Ice']=0
            tot_solid=sum(prod['Solid phases'].values())
            
            while np.abs(tot_solid) < 0.001:
                m_evap=m_evap+10
                prod =aq12.get_alldata(composition[0],composition[1],composition[2],composition[3],1000-m_evap,temp)
                if prod['Error message'] == 0:    
                    prod['Solid phases']['Ice']=0
                    tot_solid=sum(prod['Solid phases'].values())
                else:
                    tot_solid =10
            m_evap = m_evap -10   
            prod =aq12.get_alldata(composition[0],composition[1],composition[2],composition[3],1000-m_evap,temp)
            prod['Solid phases']['Ice']=0
            tot_solid=sum(prod['Solid phases'].values())
            
            
            while np.abs(tot_solid) < 0.001:
                m_evap=m_evap+1
                prod =aq12.get_alldata(composition[0],composition[1],composition[2],composition[3],1000-m_evap,temp)
                if prod['Error message'] == 0:    
                    prod['Solid phases']['Ice']=0
                    tot_solid=sum(prod['Solid phases'].values())
                else:
                    tot_solid =10
            m_evap = m_evap -1   
            prod =aq12.get_alldata(composition[0],composition[1],composition[2],composition[3],1000-m_evap,temp)
            prod['Solid phases']['Ice']=0
            tot_solid=sum(prod['Solid phases'].values())
            
            
            while np.abs(tot_solid) < 0.001:
                m_evap=m_evap+0.1
                prod =aq12.get_alldata(composition[0],composition[1],composition[2],composition[3],1000-m_evap,temp)
                if prod['Error message'] == 0:    
                    prod['Solid phases']['Ice']=0
                    tot_solid=sum(prod['Solid phases'].values())
                else:
                    tot_solid =10
            m_evap = m_evap -0.1   
            prod =aq12.get_alldata(composition[0],composition[1],composition[2],composition[3],1000-m_evap,temp)
            prod['Solid phases']['Ice']=0
            tot_solid=sum(prod['Solid phases'].values())
            
            while np.abs(tot_solid) < 0.001:
                m_evap=m_evap+0.01
                prod =aq12.get_alldata(composition[0],composition[1],composition[2],composition[3],1000-m_evap,temp)
                if prod['Error message'] == 0:    
                    prod['Solid phases']['Ice']=0
                    tot_solid=sum(prod['Solid phases'].values())
                else:
                    tot_solid =10
            m_evap = m_evap -0.01   
            prod =aq12.get_alldata(composition[0],composition[1],composition[2],composition[3],1000-m_evap,temp)
            prod['Solid phases']['Ice']=0
            tot_solid=sum(prod['Solid phases'].values())
            
            
            water_remaining = 1000- m_evap
            new_comp = 1000*np.array(composition)/water_remaining
    else:
        new_comp = np.nan
        water_remaining = np.nan
        
        
    return (new_comp)
    
def wt2mol(wtconc,substance):
    ''' accepts wtconc in fractions (0-1) and substances
    'nacl', 'mgcl', 'cacl' or 'kcl' and returns concentration in moles/kg'''
    mol_weights = {
        'nacl': 58.4428/1000,
        'mgcl': 95.2110/1000,
        'cacl': 110.9840/1000,
        'kcl': 74.5513/1000,
        'urea': 60.06/1000,
        'kfo': 84.12/1000,
        'kac': 98.14/1000,
        'cscl': 168.36/1000,
        'mgac': 142.39/1000,
        'succrose': 342.30/1000,
        }
   
    mol_conc = (1/mol_weights.get(substance))*wtconc/(1-wtconc)
    
    return(mol_conc)
    
def mol2wt(molconc,substance):
    ''' accepts wtconc in fractions (0-1) and substances
    'nacl', 'mgcl', 'cacl' or 'kcl' and returns concentration in moles/kg'''
    mol_weights = {
        'nacl': 58.4428/1000,
        'mgcl': 95.2110/1000,
        'cacl': 110.9840/1000,
        'kcl': 74.5513/1000,
        'urea': 60.06/1000,
        'kfo': 84.12/1000,
        'kac': 98.14/1000,
        'cscl': 168.36/1000,
        'mgac': 142.39/1000,
        'succrose': 342.30/1000,
        }
    wt_conc = molconc*mol_weights.get(substance)/(molconc*mol_weights.get(substance)+1)
    
    return(wt_conc)
    
def calc_real_fp(solute, concentration):
    
    if concentration > 1:
        x = concentration/100
    else:
        x = concentration
    
    liq_lines = {
        'nacl' : -543.41*x**3-5.81*x**2-59.584*x+0.0093306,  #Fitted to data from CRC handbook
        'mgcl' : -1851.3*x**3 -26.103*x**2-56.924*x -0.016935, #Fitted to data from extracted from figure in Mellinder2007
        'cacl' : -1123.3*x**3 +64.737*x**2-58.482*x +0.35097,  #Fitted to data from CRC handbook
        'kfo' : -243.25*x**3 -4.6663*x**2-47.414*x +0.099981, #Fitted to data from extracted from figure in Mellinder2007
        'kac': -340.03*x**3 -35.197*x**2-35.186*x -0.17805, #Fitted to data from extracted from figure in Mellinder2007
        'urea' : 45.526*x**3 -51.522*x**2-26.274*x -0.095706,  #Fitted to data from CRC handbook
        'succrose': -112.7*x**4 + 71.65*x**3 -28.42*x**2 -3.708*x -0.00189,  #Fitted to data from Yong and Jones (1949)
        'cscl': -106.5*x**3 + 16.58*x**2 -22.64*x + 0.0184, #Calculated from water activity data from El Guendouzi et al (2001) combined with Equation by Ge and Wang (2009)
        'kcl': -66.13*x**2-41.27*x-0.05608, # Fitted from Hall et al (1988)
        'cma': -697.7*x**3+96.71*x**2-41.96*x+0.006225, # Fitted from data read from figure in Ketcham (1991)
        }
    eut_c = {
        'nacl' : 0.233, # According to Yatsenko, O. B. and Chudotvortsev, I. G. (2002)
        'mgcl' : 0.2101, # According to Yatsenko, O. B. and Chudotvortsev, I. G. (2002)
        'cacl' : 0.305, # According to Yatsenko, O. B. and Chudotvortsev, I. G. (2002)
        'kfo' : 0.48, # According to Mellinder 2007
        'kac': 0.45, # According to Mellinder 2007
        'urea' : 0.44, # According to CRC handbook. Not specificly given, so only approximately.
        'succrose': 0.7, # Last datapoint from Yong and JOnes, not certain.
        'cscl': 0.51, # last data point in El Guendouzi et al (2001) at 51%. Fit valid to here. Eutectic, from memory, around -22?
        'kcl': 0.1955, # From Hall et al (1988)
        'cma': 0.325, # According to Ketcham (1991)
        }
    if eut_c.get(solute)<concentration:
        fp=np.nan
    else:
        fp = liq_lines.get(solute)
        
    return(fp)
    
def calc_real_liqConc(solute, temperature):
    x = np.arange(0,0.71,0.01)

    liq_lines = {
        'nacl' : -543.41*x**3-5.81*x**2-59.584*x+0.0093306,  #Fitted to data from CRC handbook
        'mgcl' : -1851.3*x**3 -26.103*x**2-56.924*x -0.016935, #Fitted to data from extracted from figure in Mellinder2007
        'cacl' : -1123.3*x**3 +64.737*x**2-58.482*x +0.35097,  #Fitted to data from CRC handbook
        'kfo' : -243.25*x**3 -4.6663*x**2-47.414*x +0.099981, #Fitted to data from extracted from figure in Mellinder2007
        'kac': -340.03*x**3 -35.197*x**2-35.186*x -0.17805, #Fitted to data from extracted from figure in Mellinder2007
        'urea' : 45.526*x**3 -51.522*x**2-26.274*x -0.095706,  #Fitted to data from CRC handbook
        'succrose': -112.7*x**4 + 71.65*x**3 -28.42*x**2 -3.708*x -0.00189,  #Fitted to data from Yong and Jones (1949)
        'cscl': -106.5*x**3 + 16.58*x**2 -22.64*x + 0.0184, #Calculated from water activity data from El Guendouzi et al (2001) combined with Equation by Ge and Wang (2009)
        'kcl': -66.13*x**2-41.27*x-0.05608, # Fitted from Hall et al (1988)
        'cma': -697.7*x**3+96.71*x**2-41.96*x+0.006225, # Fitted from data read from figure in Ketcham (1991)
        }
    eut_T = {
        'nacl' : -21.2, # According to Yatsenko, O. B. and Chudotvortsev, I. G. (2002)
        'mgcl' : -33.5, # According to Yatsenko, O. B. and Chudotvortsev, I. G. (2002)
        'cacl' : -49.5, # According to Yatsenko, O. B. and Chudotvortsev, I. G. (2002)
        'kfo' : -51, # According to Mellinder 2007
        'kac': -52, # According to Mellinder 2007
        'urea' : -17.6, # According to CRC handbook. Not specificly given, so only approximately.
        'succrose': -19.1, # Last datapoint from Yong and JOnes, not certain.
        'cscl': -23.7, # from Monnin and Dubois (1999) /Dubois(1993)-original, should find it
        'kcl': -10.69, # From Hall et al (1988)
        'cma': -27.5, # According to Ketcham (1991)
        }
    if eut_T.get(solute)>temperature:
        liq_c=np.nan
    else:
        liq_c = np.interp(-temperature, -liq_lines.get(solute), x)
        
    return(liq_c)
    
  
def calc_icemelt_cap_aw(n_nacl, n_mgcl, n_cacl, n_kcl, m_water, temp):
    ''' Accepts moles of, in given order, NaCl, MgCl, CaCl and KCl. In addition the amount of water (from 0 to inf) and the temperature is required.
'''
    m_melt = float(0)
    a_ice = aq12.get_aw(0,0,0.1,0,100000,temp)
    print(a_ice)
    a_sol = aq12.get_aw(n_nacl,n_mgcl,n_cacl,n_kcl,m_water,temp)
    print(a_sol)
    print('apa')
    
    if a_sol-a_ice<-0.005:
            while a_sol-a_ice<-0.005:
                m_melt=m_melt+10000
                a_sol = aq12.get_aw(n_nacl,n_mgcl,n_cacl,n_kcl,m_water+m_melt,temp)
            m_melt = m_melt -10000   
            a_sol = aq12.get_aw(n_nacl,n_mgcl,n_cacl,n_kcl,m_water+m_melt,temp)        
        
            while a_sol-a_ice<-0.005:
                m_melt=m_melt+1000
                a_sol = aq12.get_aw(n_nacl,n_mgcl,n_cacl,n_kcl,m_water+m_melt,temp)
            m_melt = m_melt -1000   
            a_sol = aq12.get_aw(n_nacl,n_mgcl,n_cacl,n_kcl,m_water+m_melt,temp) 
            
            while a_sol-a_ice<-0.005:
                m_melt=m_melt+100
                a_sol = aq12.get_aw(n_nacl,n_mgcl,n_cacl,n_kcl,m_water+m_melt,temp)
            m_melt = m_melt -100   
            a_sol = aq12.get_aw(n_nacl,n_mgcl,n_cacl,n_kcl,m_water+m_melt,temp) 
            
            while a_sol-a_ice<-0.005:
                m_melt=m_melt+10
                a_sol = aq12.get_aw(n_nacl,n_mgcl,n_cacl,n_kcl,m_water+m_melt,temp)
            m_melt = m_melt -10   
            a_sol = aq12.get_aw(n_nacl,n_mgcl,n_cacl,n_kcl,m_water+m_melt,temp) 
            
            while a_sol-a_ice<-0.005:
                m_melt=m_melt+1
                a_sol = aq12.get_aw(n_nacl,n_mgcl,n_cacl,n_kcl,m_water+m_melt,temp)
            m_melt = m_melt -1 
            a_sol = aq12.get_aw(n_nacl,n_mgcl,n_cacl,n_kcl,m_water+m_melt,temp) 
            
            while a_sol-a_ice<-0.005:
                m_melt=m_melt+0.10000
                a_sol = aq12.get_aw(n_nacl,n_mgcl,n_cacl,n_kcl,m_water+m_melt,temp)
            m_melt = m_melt -0.10000   
            a_sol = aq12.get_aw(n_nacl,n_mgcl,n_cacl,n_kcl,m_water+m_melt,temp) 
            
            while a_sol-a_ice<-0.005:
                m_melt=m_melt+0.010000
                a_sol = aq12.get_aw(n_nacl,n_mgcl,n_cacl,n_kcl,m_water+m_melt,temp)
            m_melt = m_melt -0.010000   
            a_sol = aq12.get_aw(n_nacl,n_mgcl,n_cacl,n_kcl,m_water+m_melt,temp) 
        
            melt_cap=m_melt 
    else:
        melt_cap=0
        
    return (melt_cap)
    
def calc_freezepoint(n_nacl, n_mgcl, n_cacl, n_kcl, m_water):
    
    temp = 0
    SI_ice = aq12.get_SI_ice(n_nacl,n_mgcl,n_cacl,n_kcl,m_water,temp)
    
    if np.abs(SI_ice-1) > 0.001:
        while np.abs(SI_ice-1) > 0.001:
            temp = temp -10
            SI_ice = aq12.get_SI_ice(n_nacl,n_mgcl,n_cacl,n_kcl,m_water,temp)
        temp = temp +10
        SI_ice = aq12.get_SI_ice(n_nacl,n_mgcl,n_cacl,n_kcl,m_water,temp)
        
        while np.abs(SI_ice-1) > 0.001:
            temp = temp -1
            SI_ice = aq12.get_SI_ice(n_nacl,n_mgcl,n_cacl,n_kcl,m_water,temp)
        temp = temp +1
        SI_ice = aq12.get_SI_ice(n_nacl,n_mgcl,n_cacl,n_kcl,m_water,temp)
        
        while np.abs(SI_ice-1) > 0.001:
            temp = temp -0.1
            SI_ice = aq12.get_SI_ice(n_nacl,n_mgcl,n_cacl,n_kcl,m_water,temp)
        temp = temp +0.1
        SI_ice = aq12.get_SI_ice(n_nacl,n_mgcl,n_cacl,n_kcl,m_water,temp)
        
        while np.abs(SI_ice-1) > 0.001:
            temp = temp -0.01
            SI_ice = aq12.get_SI_ice(n_nacl,n_mgcl,n_cacl,n_kcl,m_water,temp)
        temp = temp +0.01

    else:
        temp=0
        
    return (temp)
        
        
def calc_min_aw(n_nacl, n_mgcl, n_cacl, n_kcl,temp):
    ''' Takes amount of chemicals and a temperature. Adds 1000g of water, and removes water until none is left (ie=1). Returns water activity and the amount of remaining water '''
    
    m_water=np.arange(10000,0,-100.0)
    ie=np.zeros_like(m_water)
    AW=np.zeros_like(m_water)
    i=0
    for m in m_water:
        eqgram, fgram, SI, ENTL, ENTS, CPL, PBUB, pH, pHC, IStr, AW[i], ie[i] =aq12.get_alldata(n_nacl, n_mgcl, n_cacl, n_kcl,m,temp)
        i=i+1
    inx = np.where(ie==1)
    if inx[0].size:
        new_ref = m_water[inx[0][0]-2]
        a_ref = AW[inx[0][0]-2]
    else:
        new_ref=m_water[-2]
        a_ref=AW[-2]
#==============================================================================
#     print "first step, 100. Inital value %f, final value %f with activity %f" %(m_water[0],new_ref,a_ref)
#     plt.figure()
#     plt.plot(m_water, AW)
#==============================================================================
    
    m_water=np.arange(new_ref+100,new_ref-100,-10)
    ie=np.zeros_like(m_water)
    AW=np.zeros_like(m_water)
    i=0
    for m in m_water:
        eqgram, fgram, SI, ENTL, ENTS, CPL, PBUB, pH, pHC, IStr, AW[i], ie[i] =aq12.get_alldata(n_nacl, n_mgcl, n_cacl, n_kcl,m,temp)
        i=i+1
    inx = np.where(ie==1)
    if inx[0].size:
        new_ref = m_water[inx[0][0]-2]
        a_ref = AW[inx[0][0]-2]
    else:
        new_ref=m_water[-2]
        a_ref=AW[-2]     
#==============================================================================
#     print "second step, 10. Inital value %f, final value %f with activity %f" %(m_water[0],new_ref,a_ref)
#     plt.figure()
#     plt.plot(m_water, AW)
#==============================================================================
    
    m_water = np.arange(new_ref+10, new_ref-10,-0.1)
    ie=np.zeros_like(m_water)
    AW=np.zeros_like(m_water)
    i=0
    for m in m_water:
        eqgram, fgram, SI, ENTL, ENTS, CPL, PBUB, pH, pHC, IStr, AW[i], ie[i] =aq12.get_alldata(n_nacl, n_mgcl, n_cacl, n_kcl,m,temp)
        i=i+1
    inx = np.where(ie==1)
    if inx[0].size:
        new_ref = m_water[inx[0][0]-2]
        a_ref = AW[inx[0][0]-2]
    else:
        new_ref=m_water[-2] 
        a_ref=AW[-2]
#==============================================================================
#     print "third step, 0.10. Inital value %f, final value %f with activity %f" %(m_water[0],new_ref,a_ref)
#     plt.figure()
#     plt.plot(m_water, AW)
#==============================================================================
    
    activity = a_ref
   # rem_water = new_ref
    
    
        
    return (activity)
    
    
def calc_min_aw2(f_nacl, f_mgcl, f_cacl, f_kcl,temp):
    ''' Accepts molar fractions of each solute, and returns the minimum possible water activity for that combination'''
    
    
#==============================================================================
#     f_nacl=0
#     f_mgcl=1
#     f_cacl=0
#     f_kcl=0
#==============================================================================
    if np.sum([f_nacl, f_mgcl, f_cacl, f_kcl])!=1:
        activity_min=np.nan
        print("Error, mole fractions does not add up")
    else:
        n_nacl=f_nacl*5
        n_mgcl=f_mgcl*5
        n_cacl=f_cacl*5
        n_kcl=f_kcl*5
        
        print("Evaluating T =%d"% temp)
        
    #==============================================================================
#     mw = np.arange(10000.0,-100,-100)
#     a=np.zeros_like(mw)
#     i=0
#     for m in mw:
#         a[i]=aq12.get_aw(n_nacl,n_mgcl,n_cacl,n_kcl,m,temp)
#         i=i+1
# 
#     a[np.abs(a)<5e-2]=np.nan  
#     a_diff=np.diff(a)
#     a_diff[np.abs(a_diff)<1e-4]=0
#     plt.figure()
#     plt.plot(mw,a)
#     plt.figure()
#     plt.plot(mw[0:-1],np.diff(a))
#     
#     inx=np.where(a_diff==0)
#     inx1st=inx[0][-1]
#     
#     mw_new=mw[inx1st]
#     activity = a[inx1st]
#     print(mw_new, activity)
#     
#     
#     mw = np.arange(mw_new,-10.0,-10.0)
#     a=np.zeros_like(mw)
#     i=0
#     for m in mw:
#         a[i]=aq12.get_aw(n_nacl,n_mgcl,n_cacl,n_kcl,m,temp)
#         i=i+1
# 
#     
#     a[np.abs(a)<5e-2]=np.nan 
#     a_diff=np.diff(a)
#     a_diff[np.abs(a_diff)<1e-4]=0
#     plt.figure()
#     plt.plot(mw,a)
#     plt.figure()
#     plt.plot(mw[0:-1],np.diff(a))
#     
#     inx=np.where(a_diff==0)
#     inx1st=inx[0][-1]
#     
#     mw_new=mw[inx1st]
#     activity = a[inx1st]
#     print(mw_new, activity)
#==============================================================================
    
    act=np.array([2.0])
    mw_init=10000
    stepsize_range=[1000, 100.0, 50.0,20.0,10.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1]
    act=np.append(act,1.0)
    
    i=1
    n=0
    stepsize=stepsize_range[i]
    while np.abs(act[i]-act[i-1])>1e-4:
        mw=np.arange(mw_init,-stepsize,-stepsize)
        
        a=np.zeros_like(mw)
        j=0
        for m in mw:
            a[j]=aq12.get_aw(n_nacl,n_mgcl,n_cacl,n_kcl,m,temp)
            j=j+1
        a_ice = aq12.get_aw(1,0,0,0,10000,temp)
        a[np.abs(a)<5e-2]=np.nan
       # print a
        a[a-a_ice>-1e-3]=np.nan
        #print a
        a_diff=np.diff(a)
        #print(a_diff)
        a_diff[np.abs(a_diff)<1e-3]=0
        inx=np.where(a_diff==0)
        if len(inx[0])>0:
            #print(len(inx[0]))
            #plt.figure()
            #plt.plot(a_diff)
            inx1st=inx[0][-1]
            
            mw_init=mw[inx1st]
            act=np.append(act,a[inx1st])
            i=i+1
            stepsize=stepsize_range[i]
            n=0
        else:
            ndx=np.where(np.isfinite(a))
            mw_init=mw[ndx[0][0]]
            #print mw_init
            act=np.append(act,act[-1]+1)
            #print(mw_init, a[inx1st])
            stepsize=stepsize_range[i+1+n]
            #print stepsize
            n=n+1

    activity_min=act[-1]    
    return(activity_min)
    
def print_phases_present(n_nacl,n_mgcl,n_cacl,n_kcl,m_water,temp):
    

        
    
        
    eqgram, fgram, SI, ENTL, ENTS, CPL, PBUB, pH, pHC, IStr, AW, ie =aq12.get_alldata(n_nacl, n_mgcl, n_cacl, n_kcl,m_water,temp)
    
    eqgram[eqgram<1e-3]=0    
    fgram[fgram<1e-3]=0
    
    if eqgram[0]>0:
        no_liq_phase = 1
    else:
        no_liq_phase = 0
    
    inx_aq=np.where(eqgram>0)
    
    inx_solid=np.where(fgram>0)
    
    print("Number of species in liquid is %d"% len(inx_aq[0]))
    print("Number of solid phases is %d"% len(inx_solid[0]))
    ara=[]
    for i in np.arange(0,len(inx_aq[0]),1):
        ara.append(numbers_to_dissolved_species(inx_aq[0][i]))
    print("The dissolved species are:")
    print(ara)
    
    apa=[]
    for i in np.arange(0,len(inx_solid[0]),1):
        apa.append(numbers_to_solid_species(inx_solid[0][i]))
    print("The solid species are:")
    print(apa)
    
    return()
    
def check_possible_solids(n_nacl, n_mgcl, n_cacl, n_kcl, temp):
    
    n_nacl=1
    n_mgcl=0
    n_cacl=0
    n_kcl=0
    temp=10

    if temp > 0:
        w_phase = 68
    else:
        w_phase = 65

    if (n_nacl != 0 and n_mgcl==0 and n_cacl==0 and n_kcl==0):
        inx = [12,13,w_phase]
    elif (n_nacl != 0 and n_mgcl!=0 and n_cacl==0 and n_kcl==0):
        inx = [12,13,28,29,30,31,32,33,w_phase]
    elif (n_nacl != 0 and n_mgcl!=0 and n_cacl!=0 and n_kcl==0):
        inx = [12,13,28,29,30,31,32,33,48,49,50,51,52,54,55,w_phase]
    elif (n_nacl != 0 and n_mgcl!=0 and n_cacl!=0 and n_kcl!=0):
        inx = [12,13,20,28,29,30,31,32,33,34,48,49,50,51,52,53,54,55,w_phase]
    elif (n_nacl != 0 and n_mgcl==0 and n_cacl!=0 and n_kcl==0):
        inx = [12,13,48,49,50,51,52,w_phase]
    elif (n_nacl != 0 and n_mgcl==0 and n_cacl==0 and n_kcl!=0):
        inx = [12,13,20,w_phase]
    elif (n_nacl == 0 and n_mgcl!=0 and n_cacl==0 and n_kcl==0):
        inx = [28,29,30,31,32,33,w_phase]
    elif (n_nacl == 0 and n_mgcl!=0 and n_cacl!=0 and n_kcl==0):
        inx = [28,29,30,31,32,33,48,49,50,51,52,54,55,w_phase]
    elif (n_nacl == 0 and n_mgcl!=0 and n_cacl==0 and n_kcl!=0):
        inx = [20,28,29,30,31,32,33,34,w_phase]
    elif (n_nacl == 0 and n_mgcl!=0 and n_cacl!=0 and n_kcl!=0):
        inx = [20,28,29,30,31,32,33,34,48,49,50,51,52,53,54,55,w_phase]
    elif (n_nacl == 0 and n_mgcl==0 and n_cacl!=0 and n_kcl==0):
        inx = [48,49,50,51,52,w_phase]
    elif (n_nacl == 0 and n_mgcl==0 and n_cacl!=0 and n_kcl!=0):
        inx = [20,48,49,50,51,52,53,w_phase]
    elif (n_nacl == 0 and n_mgcl==0 and n_cacl==0 and n_kcl!=0):
        inx = [20,w_phase]
    else:
        inx=[w_phase]
    
    print(inx)    
    
    return(inx)



def aq_specID_to_dissolved_species(argument):
        switcher = {
            0: "Water",
            1: "Na+",
            2: "K+",
            3: "Mg++",
            4: "Ca++",
            5: "H+",
            6: "None1",
            7: "Cl-",
            8: "SO4--",
            9: "HSO4-",
            10: "OH-",
            11: "None2",
        }
        return switcher.get(argument, "nothing")       
            
def phase_IDnumbers_to_solid_species(argument):
        switcher = {
            0: "None ",
            1: "None ",
            2: "None ",
            3: "None ",
            4: "None ",
            5: "None ",
            6: "None ",
            7: "None ",
            8: "None ",
            9: "None ",
            10: "None ",
            11: "None ",
            12:"NaCl                    "  ,
            13:"NaCl*2H2O               "  ,
            14:"Na2SO4                  "  ,
            15:"Na2SO4*10H2O            "  ,
            16:"NaHSO4                  "  ,
            17:"NaOH                    "  ,
            18:"NaOH*H2O                "  ,
            19:"2NaOH*7H2O              "  ,
            20:"KCl                     "  ,
            21:"K2SO4                   "  ,
            22:"KHSO4                   "  ,
            23:"K2SO4*KHSO4             "  ,
            24:"KOH                     "  ,
            25:"KOH*H2O                 "  ,
            26:"KOH*2H2O                "  ,
            27:"NaK3(SO4)2              "  ,
            28:"MgCl2                   "  ,
            29:"MgCl2*2H2O              "  ,
            30:"MgCl2*4H2O              "  ,
            31:"MgCl2*6H2O              "  ,
            32:"MgCl2*8H2O              "  ,
            33:"MgCl2*12H2O             "  ,
            34:"KCl*MgCl2*6H2O          "  ,
            35:"MgSO4                   "  ,
            36:"MgSO4*H2O               "  ,
            37:"MgSO4*6H2O              "  ,
            38:"MgSO4*7H2O              "  ,
            39:"MgSO4*12H2O             "  ,
            40:"Na2SO4*MgSO4*4H2O       "  ,
            41:"3Na2SO4*MgSO4           "  ,
            42:"2Na2SO4*2MgSO4*5H2O     "  ,
            43:"KCl*MgSO4*3H2O          "  ,
            44:"K2SO4*2MgSO4            "  ,
            45:"K2SO4*MgSO4*4H2O        "  ,
            46:"K2SO4*MgSO4*6H2O        "  ,
            47:"Mg(OH)2                 "  ,
            48:"CaCl2                   "  ,
            49:"CaCl2*H2O               "  ,
            50:"CaCl2*2H2O              "  ,
            51:"CaCl2*4H2O              "  ,
            52:"CaCl2*6H2O              "  ,
            53:"KCl*CaCl2               "  ,
            54:"2MgCl2*CaCl2*12H2O      "  ,
            55:"MgCl2*2CaCl2*6H2O       "  ,
            56:"CaSO4                   "  ,
            57:"2CaSO4*H2O              "  ,
            58:"CaSO4*2H2O              "  ,
            59:"Na2SO4*CaSO4            "  ,
            60:"K2SO4*MgSO4*2CaSO4*2H2O "  ,
            61:"K2SO4*CaSO4*H2O         "  ,
            62:"K2SO4*5CaSO4*H2O        "  ,
            63:"Ca(OH)2                 "  ,
            64:"CaCl2*Ca(OH)2*H2O       "  ,
            65:"Ice                     "  ,
            66:"H2SO4                   "  ,
            67:"HCl(g)                  "  ,
            68:"H2O(g)                  "  ,
        }
        return switcher.get(argument, "nothing")
        
def getDiffusionCoeff(solute, conc):
    ''' Gives the diffusion coefficient at a given concentration (in mol/kg) at 25degC. Data from litterature, no veird corrections'''

    diff_coeff_fit = {
        'nacl' : (-0.00198*conc**4+0.01473*conc**3+0.001117*conc**2+1.466*conc+0.02569)/(conc+0.01596), # Rard and Miller 1979 or #Chang and Myerson 1985 also have data. SHould perhaps compare?
        'cacl' : (-0.08098*conc**4+0.7612*conc**3-3.208*conc**2+30.91*conc+0.2151)/(conc**3-5.77*conc**2+28.21*conc+0.1607), # c in mol/kg water, data from Rard and Miller 1979
        'kfo' : 1.67, #Lack conc dependenance data. This from CRC handbook at inf dillution
        'mgac': 0.92, #Lack conc dependenance data. This from CRC handbook at inf dillution
        'urea' : 1.40171 - 5.91362e-2*conc - 3.90011e-4*conc**2, # From Sorell and Myerson 1982
        'succrose': (-0.0802*conc+1.148)/(conc+2.184), #English and Dole (1950) above 60wt% and Henrion (1964) below that
        'cscl': (-0.267*conc**4+5.919*conc**3+27.21*conc**2+47.1*conc+0.7114)/(conc**3+15.76*conc**2+25.37*conc+0.3489),   #Rard and Miller 1981
        'kcl': (-0.002328*conc**4 + 0.01335*conc**3 +0.07635*conc**2 + 1.799*conc +0.06447)/(conc+0.0324),   #Gosting(1950)
        'mgcl': (-2.922*conc**2+31.6*conc+0.2495)/(conc**3-6.042*conc**2+31.13*conc+0.201), #Miller et al 1984
        }
    max_C = {
        'nacl' : 6.02, #mol/kg # Rard and Miller 1979 
        'cacl' : 7, # # c in mol/kg water, data from Rard and Miller 1979
        'kfo' : 100, # 
        'mgac': 100, # 
        'urea' : 10.5, # According to CRC handbook. Not specificly given, so only approximately.
        'succrose': 8.5, ##English and Dole (1950)
        'cscl': 11.1, #  #Rard and Miller 1981
        'kcl': 4.5, # #Gosting(1950)
        'mgcl': 5.8, #Miller et al 1984
        }
    if conc > max_C.get(solute):
        diff=np.nan
    else:
        diff = diff_coeff_fit.get(solute)
        
    return(diff)

#def getSolutionDensity(solute,conc):
     #''' Gives the solution density at a given concentration at 25degC. Data from litterature, no veird corrections'''
#==============================================================================
# def plot_phasediagram(f_nacl, f_mgcl, f_cacl, f_kcl, min_temp):
#     f_nacl=1
#     f_mgcl=0
#     f_cacl=0
#     f_kcl=0
#     min_temp = -25    
#     
#     temps = np.arange(0,-25.1,-0.1)
#     conc = np.arange(0,15,0.1)
#     
#     for T in temps:
#         
#         for c in conc
#     
#     
#     eqgram, fgram, SI, ENTL, ENTS, CPL, PBUB, pH, pHC, IStr, AW, ie =aq12.get_alldata(n_nacl, n_mgcl, n_cacl, n_kcl,m_water,temp)
#==============================================================================
def calc_MC_from_phase(self, temp, conc, msol, solute):
     cf= calc_real_liqConc(solute,temp)
     MC = (conc*msol)/cf - msol
     return(MC)




def find_precip_temp(n_nacl, n_mgcl, n_cacl, n_kcl,m_water):
    ''' Composition is 4x1 vectors giving the molar concentration of
NaCl, MgCl, CaCl, KCl in the given order. Temp begins at 200 degC. The script begins with the given composition dissolved in a 1000 gram of water. 
If the given composition does not fully dissolve in 1000g of water, the function returns nan. This should be adjusted when I come up with a good idea of what is most relevant '''
    
#    n_nacl=5
#    n_mgcl=3
#    n_cacl=4
#    n_kcl=0
#    m_water=2000
    def increment_property(temp_initial,increment, prev_dict):
        prod =aq12.get_alldata(composition[0],composition[1],composition[2],composition[3],m_water,temp_initial)
        d=dict((k,temp_initial) for k, v in prod['Solid phases'].items() if v > 0)
        #print('The initial temperature for the loop in %.2f' %temp_initial)
        #print('And the increment is %.2f' %increment)
        #print('The sum of the two are %.2f' %(temp_initial+increment))
        #print('The initial dict keys for the loop are')
        print(d.keys())
        while (d.keys() == prev_dict.keys() and prod['Error message'] == 0):
            #print('The step temp should be %.2f'% (temp_initial+increment))
            temp_initial=(temp_initial+increment)
            #print('but step temp is %.2f'% temp_initial)
            prod =aq12.get_alldata(composition[0],composition[1],composition[2],composition[3],m_water,temp_initial)
            d=dict((k,temp_initial) for k, v in prod['Solid phases'].items() if v > 0)
            #print('The new d:')
            #print(d.keys())
            #print('At temperarure %.2f' %temp_initial)
            #print('Which is compared with:')
            #print(prev_dict.keys())
            #print('-----------------')
        return(temp_initial-increment,d, prod['Error message'])
        
    composition = [n_nacl, n_mgcl, n_cacl, n_kcl]
    T_init=50.05
    data_dict={}
    prod =aq12.get_alldata(composition[0],composition[1],composition[2],composition[3],m_water,T_init)
    if prod['Error message'] == 0:
        d=dict((k,T_init) for k, v in prod['Solid phases'].items() if v > 0)
        for k,v in d.items():
            if k not in data_dict: 
                data_dict[k]=v

        prev_dict = data_dict
        prd_error=0

        while prd_error == 0:
            print('.....................................................')
            print('Start loop from temperature %.2f and error is %f' %(T_init-0.05, prd_error))
            print('With new reference dict being:')
            print(prev_dict.keys())
            T_init,di,prd_error=increment_property(T_init-0.05,-10,prev_dict)
            #print(ara)
            T_init,di,prd_error=increment_property(T_init,-5,prev_dict)
            #print(ara)
            T_init,di,prd_error=increment_property(T_init,-1,prev_dict)
            #print(ara)
            T_init,di,prd_error=increment_property(T_init,-0.05,prev_dict)
            for k,v in di.items():
                if k not in data_dict: 
                    data_dict[k]=v
            prev_dict = data_dict
         
    return (data_dict)