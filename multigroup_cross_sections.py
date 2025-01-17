# -*- coding: utf-8 -*-
"""
Multigroup Photon Cross Sections

Autores: E. Marqués
Fecha:   11/01/2025
"""

import numpy as np
import pandas as pd

class MultigroupPhotonCrossSections:
    
    def __init__(self):
        self.electron_classical_radius_squared_barn = 0.07939827
        self.electron_rest_mass_keV = 511.006
        self.group_start = (150,125,100,75,50,25)
        self.group_stop = (125,100,75,50,25,0)
        self.group_step = (-1,-1,-1,-1,-1,-1)
# Biggs F, Lighthill R. Analytical approximations for x-ray cross sections III. 
# Sandia Natl. Lab., vol. SAND87, no. 70; 1988.        
        
        self.element_data = {
            'H':{'Z':1,'A':1.00794,'Conv':1.67,'Z/A':0.9921},\
            'O':{'Z':8,'A':15.9994,'Conv':26.57,'Z/A':0.5000},\
            'Al':{'Z':13,'A':26.98154,'Conv':44.80,'Z/A':0.4818}
            }
        #self.compound_data = {
        #    'Water':{'H':2.*self.element_data['H']['A']/\
        #             (2.*self.element_data['H']['A']+self.element_data['O']['A']),\
        #             'O':self.element_data['O']['A']/\
        #             (2.*self.element_data['H']['A']+self.element_data['O']['A'])}
        #    }
        self.photoelectric_cross_sections_fitting_parameters = {
            'H':pd.DataFrame(
                data=[[0.01,0.014,1.000E-08,0.,0.,0.],\
                      [0.014,0.1,-6.383E+01,-6.446E+00,1.317E+01,-5.045E-02],\
                      [0.1,0.8,3.051E+00,-7.818E+00,1.144E+01,6.959E-02],\
                      [0.8,4.,7.636E-02,-9.406E-01,6.144E+00,1.425E+00],\
                      [4.,20.,1.180E-03,-8.236E-02,2.886E+00,5.534E+00],\
                      [20.,100.,1.620E-05,-5.610E-03,1.214E+00,1.761E+01],\
                      [100.,500.,1.034E-06,-4.114E-04,6.287E-01,3.927E+01]],\
                index=[1,2,3,4,5,6,7],\
                columns=['START','FINISH','A_1','A_2','A_3','A_4']),\
            'O':pd.DataFrame(
                data=[[0.01,0.0483,1.144E+04,0.,0.,0.],\
                      [0.0483,0.532,-2.863E+02,4.085E+02,4.436E+01,-1.782E+00],\
                      [0.532,4.,-7.181E+01,4.748E+02,5.542E+03,-1.363E+03],\
                      [4.,20.,2.745E+00,-1.747E+02,7.159E+03,-2.213E+03],\
                      [20.,100.,3.774E-02,-1.559E+01,4.045E+03,1.810E+04],\
                      [100.,500.,3.169E-03,1.473E+00,7.214E+02,4.048E+05]],\
                index=[1,2,3,4,5,6],\
                columns=['START','FINISH','A_1','A_2','A_3','A_4']),\
            'Al':pd.DataFrame(
                data=[[0.01,0.0159,-1.654E+04,1.585E+02,3.907E+00,-3.383E-02],\
                      [0.0159,0.073,1.122E+03,-4.015E+01,6.623E-01,-2.813E-03],\
                      [0.073,0.1177,2.390E+04,-6.953E+02,-7.978E+01,1.974E+00],\
                      [0.1177,1.560,-5.284E+02,1.399E+03,4.360E+02,-4.747E+01],\
                      [1.560,20.,-3.674E+00,-1.622E+01,2.732E+04,-1.752E+04],\
                      [20.,100.,4.158E-01,-1.351E+02,2.716E+04,3.723E+02],\
                      [100.,500.,1.125E-02,2.747E+00,1.174E+04,5.695E+05]],\
                index=[1,2,3,4,5,6,7],\
                columns=['START','FINISH','A_1','A_2','A_3','A_4']) 
                }

# return the mass absorption coefficient of an element [cm^2 g^-1] 
# element: chemical symbol of the element ['H','O']
# E:       energy bin [keV] 
    def get_mass_absorption_coefficient_element(self,element,E):
        df = self.photoelectric_cross_sections_fitting_parameters[element]
        A_1 = df[(df['START']<=E) & (df['FINISH']>= E)]['A_1'].values[0]
        A_2 = df[(df['START']<=E) & (df['FINISH']>= E)]['A_2'].values[0]
        A_3 = df[(df['START']<=E) & (df['FINISH']>= E)]['A_3'].values[0]
        A_4 = df[(df['START']<=E) & (df['FINISH']>= E)]['A_4'].values[0]
        mass_absorption_coefficient = A_1/E+A_2/E**2+A_3/E**3+A_4/E**4
        return mass_absorption_coefficient

# return the mass scattering coefficient of an element [cm^2 g^-1] 
# element: chemical symbol of the element ['H','O']
# E:       energy bin [keV]    
    def get_mass_scattering_coefficient_element(self,element,E):
        df = self.element_data[element]
        Z = df['Z']
        A = df['A'] # g mole^-1
        L = 0.40061 # cm^2 mole^-1
        X = E/511.04
        R = L*(1.+1.148*X+0.06141*X**2)/(1.+3.171*X+0.9328*X**2+0.02572*X**3)
        mass_scattering_coefficient = R*(Z/A)
        return mass_scattering_coefficient
    
# return the Klein-Nishina mass scattering coefficient of an element [cm^2 g^-1] 
# element: chemical symbol of the element ['H','O']
# E:       energy bin [keV]    
    def get_KN_mass_scattering_coefficient_element(self,element,E):
        df = self.element_data[element]
        Z = df['Z']
        A = df['A'] # g mole^-1
        L = 0.40061 # cm^2 mole^-1
        X = E/511.04
        omega = 1./(1.+2*X)
        R = (3/4)*L*(((2.+2.*X-X**2)/(2.*X**3))*np.log(omega)+\
                     2.*omega*(1.+X)**2/X**2-omega**2*(1.+3*X))
        mass_scattering_coefficient = R*(Z/A)
        return mass_scattering_coefficient
    
# return the mass attenuation coefficient of an element [cm^2 g^-1] 
# element: chemical symbol of the element ['H','O']
# E:       energy bin [keV]    
    def get_mass_attenuation_coefficient_element(self,element,E):    
        mass_absorption_coefficient = \
            self.get_mass_absorption_coefficient_element(element,E)
        mass_scattering_coefficient = \
            self.get_mass_scattering_coefficient_element(element,E)
        mass_attenuation_coefficient = \
            mass_absorption_coefficient+mass_scattering_coefficient
        return mass_attenuation_coefficient
    
# return the mass absorption coefficient of a compound [cm^2 g^-1] 
# compound: weight fraction of each element in the compound
# E:       energy bin [keV] 
    def get_mass_absorption_coefficient_compound(self,compound,E):
        pass
        
        
    def get_group_angle_transfer_matrix_element(self,group_in,group_out,direction_in,direction_out):        
        chi_m = np.dot(direction_in,direction_out)
        group_angle_transfer_matrix_element = 0.
        group_in_start = self.group_start[group_in]
        group_in_stop = self.group_stop[group_in]
        group_in_step = self.group_step[group_in]
        group_out_start = self.group_start[group_out]
        group_out_stop = self.group_stop[group_out]
        group_out_step = self.group_step[group_out]
        integral_in = 0.
        for E_in in np.arange(group_in_start,group_in_stop,group_in_step):
            #print('E_in=',E_in)
            lambda_in = self.electron_rest_mass_keV/E_in
            integral_out = 0.
            for E_out in np.arange(group_out_start,group_out_stop,group_out_step):
                #print('E_out=',E_out)
                E_max = E_in
                E_min = E_in/(1.+2./lambda_in)
                #print('E_max=',E_max)
                #print('E_min=',E_min)
                if((E_out<=E_max)&(E_out>=E_min)):
                    print('E_in=',E_in,'E_out=',E_out,'E_max=',E_max,'E_min=',E_min)
                    lambda_out = self.electron_rest_mass_keV/E_out
                    chi = 1.+lambda_in-lambda_out
                    print('chi=',chi)
                    if(chi==chi_m):
                        integral_out += E_in/E_out+E_out/E_in-1.+chi**2
            integral_in += group_in_step/(E_in**2*(group_in_start-group_in_stop))*integral_out
        group_angle_transfer_matrix_element = 0.5*self.electron_classical_radius_squared_barn*self.electron_rest_mass_keV*integral_in    
        return group_angle_transfer_matrix_element
        
        
        
        
        