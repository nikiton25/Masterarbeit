# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 11:27:55 2022

@author: nt
"""
import math 
import numpy as np
import pandas as pd
# from thermo.chemical import Chemical
# from thermo.chemical import Mixture
# from scipy.interpolate import RegularGridInterpolator
from tabulate import tabulate


# convert into dataframe
df = pd.read_excel(r'C:\Users\nt\modell-waermeuebertrager\modell-waermeuebertrager\Eingabemaske_Plattenwärmeübertrager.xlsx', sheet_name='Plattenwärmeübertrager', index_col=0, header=0, usecols="B:C")
df = df.dropna()
dict = df.to_dict()
print(dict)

print(dict.keys())
print(dict.values())

T1 = float(df['Value']['T1'])
T2 = float(df['Value']['T2'])
t1 = float(df['Value']['t1'])
t2 = float(df['Value']['t2'])

class Plateheatexchanger():

    def __init__(self, fc, m, beta, Nt, Np, s, Lw, Lv, Lh, Lc, Dp, k_steel,  PD_Design_prim, PD_Design_sec, Q, A_tot):
                
        self.flow_condition          =   fc    #[-] Excel data
        self.Medium                  =   m     #[-] Excel data
        self.Chevron_angle           =   beta     #[°] Excel Data
        self.Tot_No_Plates           =   Nt     #[-] Excel Data
        self.No_of_passes            =   Np     #[-] Excel Data

        self.Plate_Thickness         =   s     #[m] Excel Data
        self.Horozontal_plate_width  =   Lw     #[m] Excel Data
        self.Vertical_port_length    =   Lv     #[m] Excel Data
        self.Horizontal_port_width   =   Lh     #[m] Excel Data
        self.Compressed_plate_length =   Lc     #[m] Excel Data
        self.Port_diameter           =   Dp     #[m] Excel Data
 
        self.Plate_thermal_cond_SS   =   k_steel     #[W/m,°C] Excel Data

        self.Pressure_Drop_Prim      =   PD_Design_prim     #[kPa] Excel Data
        self.Pressure_Drop_Sec       =   PD_Design_sec     #[kPa] Excel Data

        self.Heat_load               =   Q  #[W] user defined
        self.A_tot_eff               =   A_tot  #[m2] Datenblatt
        
        # 'Duty Requirements'
        # 'Primary circuit'
        # self.Inlet_Temp_prim         =   T1          #[K] Excel Data
        # self.Outlet_Temp_prim        =   T2       #[K] Excel Data
        # 'Secondary circuit'
        # self.Inlet_Temp_sec          =   t1       #[K] Excel Data
        # self.Outlet_Temp_sec         =   t2         #[K] Excel Data
        # 'Reference Temperature'
        # self.T_ref_h                 =  (T1 + T2) / 2
        # self.T_ref_c                 =  (t1 + t2) / 2

        
        # self.calculator()
        
    
    def calculator(self):
        global T1, T2, t1, t2
        
        self.Fouling(),
        self.TemperatureCalculations(T1, T2, t1, t2),
        self.PhysicalProperties(),
        self.physicalPlateCalculations(),
        self.flowRelatedCalculation(T1, T2, t1, t2),
        self.ReserveCalculation(),
        self.pressureLossCalculation(),
        self.NTUMethod(T1, T2, t1, t2),
        self.print_output(T1, T2, t1, t2)
        self.getNewTemp(T1, T2, t1, t2)
        return

    def Fouling(self):
        if self.Medium =='water':
           self.Rdi = 0.00002 #[m²,°C/W]
           self.Rdo = 0.00002 #[m²,°C/W]
        elif self.Medium =='gas':
           self.Rdi = 0.00002 #[m²,°C/W]
           self.Rdo = 0.00005 #[m²,°C/W]
        else:
           self.Rdi = 0.00002 #[m²,°C/W]
           self.Rdo = 0.00008 #[m²,°C/W]
        return self.Rdi, self.Rdo
    
    def TemperatureCalculations(self, T1, T2, t1, t2):
        self.T_ref_h                 =  (T1 + T2) / 2
        self.T_ref_c                 =  (t1 + t2) / 2
        dT1 = T1-t2 #[K]
        dT2 = T2 - t1 #[K]
        self.dTlm = (dT1-dT2)/math.log(dT1/dT2)
        return self.dTlm, self.T_ref_h, self.T_ref_c
    
    def PhysicalProperties(self):
        'Primary circuit'
        self.Den_water_prim          =   966.2  #[kg/m³] user defined
        self.Kin_Vis_water_prim    =   ((2.414*10**-5)*10**(247.8/(self.T_ref_h-142))/self.Den_water_prim) #[m2/s] Excel Data
        # self.Kin_Vis_water_prim      =   0.330159*10**-6 #[m2/s] Excel Data
        self.Kin_Vis_water_wall_prim =   0.4067*10**-6   #[m2/s] user defined    
        self.thermal_cond_water_prim = ((-8.354 * 10**-6)*(self.T_ref_h**2) + (6.53 * 10**-3) * (self.T_ref_h) - 0.5981)  #[W/m,°C] 
        # self.thermal_cond_water_prim = 0.6746
        self.Heat_cap_water_prim     =   4206         #[J/kg,°C] user defined
        '''Secondary circuit'''
        self.Den_water_sec           =   981.7  #[kg/m³] user defined  
        self.Kin_Vis_water_sec       =   ((2.414*10**-5)*10**(247.8/(self.T_ref_c-142))/self.Den_water_sec) #[m2/s] Excel Data
        # self.Kin_Vis_water_sec       =   0.45635*10**-6 #[m2/s] Excel Data
        self.Kin_Vis_water_wall_sec  =   0.4105*10**-6  #[m2/s] user defined    
        self.thermal_cond_water_sec  =   ((-8.354 * 10**-6)*(self.T_ref_c**2) + (6.53 * 10**-3) * (self.T_ref_c) - 0.5981)  #[W/m,°C] 
        # self.thermal_cond_water_sec = 0.6570
        self.Heat_cap_water_sec      =   4187     #[J/kg,°C] user defined
        return
    
    def physicalPlateCalculations(self):
        self.Ne = self.Tot_No_Plates - 2
        self.p = (self.Compressed_plate_length/self.Tot_No_Plates)
        self.b = self.p - self.Plate_Thickness
        self.A_ch = self.b * self.Horozontal_plate_width
        self.A_1 = self.A_tot_eff / self.Ne
        self.Lp = self.Vertical_port_length - self.Port_diameter
        self.A_1p = self.Lp * self.Horozontal_plate_width
        self.enl_factor = self.A_1 / self.A_1p
        self.Dh = (2 * self.b) / self.enl_factor
        self.N_cp = (self.Tot_No_Plates - 1) / (2* self.No_of_passes)
        self.A_port = ((math.pi/4) * self.Port_diameter**2)
        return self.Ne, self.p, self.b, self.A_ch, self.A_1, self.Lp,self.A_1p, self.enl_factor,self.Dh,  self.N_cp, self.A_port

    def flowRelatedCalculation(self, T1, T2, t1, t2):
        self.m_prim = self.Heat_load/(self.Heat_cap_water_prim*(T1 - T2))
        self.m_sec = self.Heat_load/(self.Heat_cap_water_sec*(t2 - t1))
        self.u_ch_h = self.m_prim / (self.Den_water_prim*self.A_ch*self.N_cp)
        self.u_ch_c = self.m_sec / (self.Den_water_sec*self.A_ch*self.N_cp)
        self.G_ch_h = self.m_prim / (self.N_cp * self.A_ch)
        self.G_ch_c = self.m_sec / (self.N_cp * self.A_ch)
        self.G_port_h = self.m_prim /self.A_port
        self.G_port_c = self.m_sec / self.A_port
        self.u_p_h = self.m_prim / (self.Den_water_prim*self.A_port)
        self.u_p_c = self.m_sec / (self.Den_water_sec*self.A_port)
        return self.m_prim, self.m_sec, self.u_ch_h, self.u_ch_c, self.G_ch_h, self.G_ch_c,self.G_port_h, self.G_port_c,self.u_p_h,  self.u_p_c

    def ReserveCalculation (self):
        self.Pr_h = (self.Heat_cap_water_prim * self.Kin_Vis_water_prim * self.Den_water_prim) / self.thermal_cond_water_prim
        self.Pr_c = (self.Heat_cap_water_sec * self.Kin_Vis_water_sec * self.Den_water_sec) / self.thermal_cond_water_sec
        self.Re_h = (self.u_ch_h * self.Dh) / self.Kin_Vis_water_prim
        self.Re_c = (self.u_ch_c * self.Dh) / self.Kin_Vis_water_sec
        self.Nu_h = 0.36* (self.Re_h**(2/3)) * (self.Pr_h**(1/3))
        self.Nu_c = 0.36* (self.Re_c**(2/3)) * (self.Pr_c**(1/3))
        self.alpha_h=(self.Nu_h*self.thermal_cond_water_prim)/ self.Dh
        self.alpha_c=(self.Nu_c*self.thermal_cond_water_sec)/ self.Dh
        
        self.Uo = 1.0/((1.0/self.alpha_h)+(1.0/self.alpha_c)+(self.Plate_Thickness/self.Plate_thermal_cond_SS))
        self.Uf= 1.0/((1.0/self.alpha_h)+(1.0/self.alpha_c)+(self.Plate_Thickness/self.Plate_thermal_cond_SS)+ self.Rdi + self.Rdo)
        
        self.Qc = self.Uo * self.A_tot_eff * self.dTlm
        self.Qf = self.Uf * self.A_tot_eff * self.dTlm 
        
        self.Res_Q = ((self.Qc/self.Heat_load)-1)*100 
        return self.Pr_h, self.Pr_c, self.Re_h, self.Re_c, self.Nu_h, self.Nu_c,self.alpha_h, self.alpha_c,self.Uo,  self.Uf, self.Qc, self.Qf, self.Res_Q

    def FrictionFactor(self): #Kumar Korelation
        if self.Re_h < 100:
            self.f_h = 19.400 / (self.Re_h**0.589)
            self.f_c = 19.400 / (self.Re_c**0.589)
        else:
            self.f_h = 2.99 / (self.Re_h**0.183)
            self.f_c = 2.99 / (self.Re_c**0.183)
        return self.f_h, self.f_c  
    
    def pressureLossCalculation(self):
        self.P_p_h = 1.4 * self.No_of_passes * ((self.G_port_h**2)/(2*self.Den_water_prim))
        self.P_p_c = 1.4 * self.No_of_passes * ((self.G_port_c**2)/(2*self.Den_water_sec))
        self.FrictionFactor()
        self.P_c_h = 4 * self.f_h * ((self.Vertical_port_length * self.No_of_passes) / self.Dh) * (self.G_ch_h **2 / (2 * self.Den_water_prim)) * (self.Kin_Vis_water_wall_prim / self.Kin_Vis_water_prim)**-0.17 
        self.P_c_c = 4 * self.f_c * ((self.Vertical_port_length * self.No_of_passes) / self.Dh) * (self.G_ch_c **2 / (2 * self.Den_water_sec)) * (self.Kin_Vis_water_wall_sec / self.Kin_Vis_water_sec)**-0.17
        self.P_t_h = self.P_p_h + self.P_c_h
        self.P_t_c = self.P_p_c + self.P_c_c
        return 
    
    def NTUMethod(self, T1, T2, t1, t2):
        self.C_h = self.Heat_cap_water_prim*self.m_prim
        self.C_c = self.Heat_cap_water_sec*self.m_sec
        self.C_min = min(self.C_h,self. C_c)
        self.C_max = max(self.C_h, self.C_c)
        self.NTU = (self.Uo*self.A_tot_eff) / (self.C_min)
        self.R = self.C_min / self.C_max
        if self.flow_condition =='counterflow':
            self.P_cc = (1-(np.exp(-self.NTU*(1-self.R)))) / (1-(self.R*(np.exp(-self.NTU*(1-self.R))))) #for R:=1 counterflow
            self.Q_new = self.P_cc * self.C_min * (T1 - t1) 
        else:
            self.P_co = (1-np.exp(-self.NTU(1+self.R))) / (1+self.R) #for R:=1 co-current
            self.Q_new = self.P_co * self.C_min * (T1 - t1) 
        
        self.T2_new = T1 - (self.Q_new / self.C_h) 
        self.t2_new = t1 + (self.Q_new / self.C_c)
        return

    def print_output(self, T1, T2, t1, t2):
        
        print(tabulate([['Flow Condition', self.flow_condition],
                        ['Medium', self.Medium ], 
                        ['Chevron angle', self.Chevron_angle], 
                        ['Total No of Plates', self.Tot_No_Plates], 
                        ['Total No of Passes', self.No_of_passes],
                        ['Plate_Thickness', self.Plate_Thickness],
                        ['Horozontal_plate_width', self.Horozontal_plate_width],
                        ['Vertical_port_length', self.Vertical_port_length],
                        ['Horizontal_port_width', self.Horizontal_port_width],   
                        ['Compressed_plate_length', self.Compressed_plate_length],
                        ['Port_diameter', self.Port_diameter],
                        ['Plate_thermal_cond_SS', self.Plate_thermal_cond_SS],
                        ['Pressure_Drop_Prim', self.Pressure_Drop_Prim],   
                        ['Pressure_Drop_Sec', self.Pressure_Drop_Sec],
                        ['Heat_load', self.Heat_load],
                        ['A_tot_eff', self.A_tot_eff],
                        ['Inlet_Temp_prim', T1],   
                        ['Outlet_Temp_prim', T2],  
                        ['Inlet_Temp_sec', t1],
                        ['Outlet_Temp_sec', t2],
                        ['Referenz Temp Prim', self.T_ref_h],
                        ['Referenz Temp Sec', self.T_ref_c],
                        ['Pitch', self.p], 
                        ['mean channel spacing', self.b], 
                        ['channel flow Area', self.A_ch], 
                        ['enlargment factor', self.enl_factor], 
                        ['hydraulic diameter', self.Dh],
                        ['Inner Fouling', self.Rdi],
                        ['Outer Fouling', self.Rdo],
                        ['flow channel velocity hot streams', self.u_ch_h], 
                        ['flow channel velocity cold streams', self.u_ch_c],
                        ['mass velocity per channel hot streams', self.G_ch_h],
                        ['mass velocity per channel cold streams', self.G_ch_c],
                        ['port area', self.A_port],
                        ['mass velocity per port hot streams', self.G_port_h],
                        ['mass velocity per port cold streams', self.G_port_c],
                        ['port channel velocity hot streams', self.u_p_h],
                        ['port channel velocity hot streams', self.u_p_c],
                        ['Prandl No hot streams', self.Pr_h],
                        ['Prandl No cold streams', self.Pr_c],
                        ['Reynolds No hot streams', self.Re_h],
                        ['Reynolds No cold streams', self.Re_c],
                        ['Nusselt No hot streams', self.Nu_h],
                        ['Nusselt No cold streams', self.Nu_c],
                        ['heat transfer coeff hot streams', self.alpha_h],
                        ['heat transfer coeff cold streams', self.alpha_c],
                        ['clean heat transfer coeff', self.Uo],
                        ['fouled heat transfer coeff', self.Uf],
                        ['heat duties cleaned', self.Qc],
                        ['heat duties fouled', self.Qf],
                        ['Reserve', self.Res_Q],
                        ['pressure drop port hot stream', self.P_p_h],
                        ['pressure drop channel hot stream', self.P_c_h],
                        ['total pressure drop hot stream', self.P_t_h],
                        ['pressure drop port cold stream', self.P_p_c],
                        ['pressure drop channel cold stream', self.P_c_c],
                        ['total pressure drop cold stream', self.P_t_c],
                        ['NTU', self.NTU],
                        ['Verhältnis der Wärmekapazitätsströme R', self.R],
                        ['Heat duty new', self.Q_new],
                        ['Prim RL Temp', self.T2_new],
                        ['Sec RL Temp', self.t2_new],
                        # ['Min VL Temp', self.T1_min]
                        ], headers=['Paramter', 'Wert']))
        
        
        

            

    # def getNewTemp(self, T1, T2, t1, t2):
        
    #     
    #     while (self.Res_Q >= 0):
                # T1 = T1 - 1
    #           print('The minimum inlet Temperature for this Heat exchanger is ', T1, 'K')            
    #       else:
    #         print ('Increase your Heat Exchanger')
    #     return T1
            
    # def getNewTemp(self, T1, T2, t1, t2):
        
    #     self.Res_Q >= 0
    #     while True:
    #         try:
    #             T1 = T1 - 1
    #             print('The minimum inlet Temperature for this Heat exchanger is ', T1, 'K')
    #         except StopIteration:
    #             self.Res_Q < 0
    #             break
                
    # def getNewTemp(self, T1, T2, t1, t2):
        
    #     if self.Res_Q >= 0:
    #         for i in range (1):
    #             T1 = T1 - i
    #             print('The minimum inlet Temperature for this Heat exchanger is ', T1, 'K')
    #         else:
    #             print ('Increase your Heat Exchanger')
    #     return T1
            
    # def getNewTemp(self, T1, T2, t1, t2):
        
    #     if self.Res_Q >= 0:
    #         while (T1 >= T2):
    #             T1 = T1 - 1
    #             print('The minimum inlet Temperature for this Heat exchanger is ', T1, 'K')
    #     else:
    #             print ('Increase your Heat Exchanger')
    #     return T1    
    
    # def getNewTemp(self, T1, T2, t1, t2):
    #     while (self.Res_Q >= 0):
    #         T1-= 1 
    #         if T1 >= T2:
    #             print('The minimum inlet Temperature for this Heat exchanger is ', T1, 'K')
    #             break
    #         else:
    #             print ('Increase your Heat Exchanger')
    #     return T1        
    
    # for T1 in enumerate(T1)
    
    # def getNewTemp(self, T1, T2, t1, t2):
    #     for T1 in self.calculator():
    #         T1-= 1 
    #         if self.Res_Q >= 0 :
    #             print('The minimum inlet Temperature for this Heat exchanger is ', T1, 'K')
    #             continue
    #         if self.Res_Q <= 0 :
    #             break
    #         print ('Increase your Heat Exchanger')
    #         # else:
    #         #     print ('Increase your Heat Exchanger')
   
    def getNewTemp(self, T1, T2, t1, t2):
        while self.Res_Q >= 0:
            if T1<T2:
                T1 -=1
                print('The minimum inlet Temperature for this Heat exchanger is ', T1, 'K')
            else:
                print ('Increase your Heat Exchanger')
                



    
    # def get_output(self):
    #         resluts_dict = {
    #         "Reserve Heatexchanger"    : self.Res_Q
    #         }
    #         return resluts_dict
                
  
newClass = Plateheatexchanger(df['Value']['fc'], 
                              df['Value']['m'], 
                              float(df['Value']['beta']), 
                              float(df['Value']['Nt']), 
                              float(df['Value']['Np']), 
                              float(df['Value']['s']), 
                              float(df['Value']['Lw']),
                              float(df['Value']['Lv']),
                              float(df['Value']['Lh']),
                              float(df['Value']['Lc']), 
                              float(df['Value']['Dp']), 
                              float(df['Value']['k_steel']), 
                              float(df['Value']['PD_Design_prim']),
                              float(df['Value']['PD_Design_sec']),
                              float(df['Value']['Q']),
                              float(df['Value']['A_tot']))







newClass.calculator()


# newClass = Plateheatexchanger(fc, m, beta, Nt, Np, s, Lw, Lv, Lh, Lc, Dp, k_steel,  PD_Design_prim, PD_Design_sec, Q, 
#              A_tot, T1, T2, t1, t2)
# newClass = Plateheatexchanger(fc=df[], m, beta, Nt, Np, s, Lw, Lv, Lh, Lc, Dp, k_steel,  PD_Design_prim, PD_Design_sec, Q, 
#              A_tot, T1, T2, t1, t2)

# val = newClass.ReserveCalculation()
# print(val)



# Parameters = list(df.Parameters)
# Symbol = list(df.Symbol)
# Value = list(df.Value)
# Dimensions = list(df.Dimensions)


# Werte = []

# for Parameters in Parameters:
#     for Symbol in Symbol:
#         for Value in Value: 
#             for Dimensions in Dimensions:
#                 Werte.append([Parameters, Symbol, Value, Dimensions])
                
# print(Werte)

# def main ():
#     # HX_1021 = Plateheatexchanger(fc='Counterflow', m=df['Value]['m'], 45, 22, 1, 0.006, 1.19, 2.43, 0.72, 0.533, 0.24, 16.5,15, 15, 43000, 0.62, 120, 57.50, 55, 70.62)

# #     print(HX_1021)
    
        
        
# df['Value']['m']      
# print        
        
        
        
        
        
        
        
        
        # output = {
        #     "Fouling"  : self.Fouling(),
        #     "MassFlow"  : self.TemperatureCalculations(),
        #     "dTLm"     : self.physicalPlateCalculations(),
        #     "Eff_No_Plates"       : self.flowRelatedCalculation(),
        #     "Pitch"        : self.ReserveCalculation()
        #     # "Mean_channel_spacing"       : self.pressureLossCalculation(),
        #     # "Channel_flow_area"       : self.Channel_flow_area(),     
        #     # "Single_plate_heat_transfer_area"  : self.Single_plate_heat_transfer_area(),
        #     # "Projected_plate_area"  : self.Projected_plate_area(),
        #     # "Enlargment_factor"     : self.Enlargment_factor(),
        #     # "Hydraulic_diameter_flow_channel"       : self.Hydraulic_diameter_flow_channel(),
        #     # "Number_channels_per_pass"        : self.Number_channels_per_pass(),
        #     # "Channel_velocity"       : self.Channel_velocity(),
        #     # "Mass_velocity_per_channel"       : self.Mass_velocity_per_channel,
        #     # "Mass_velocity_port"  : self.Mass_velocity_port(),
        #     # "Port_velocity"  : self.Port_velocity(),
        #     # "Prandl"     : self.Prandl(),
        #     # "Reynolds"       : self.Reynolds(),
        #     # "Nusselt"        : self.Nusselt(),
        #     # "Heat_transfer_coeff"       : self.Heat_transfer_coeff(),
        #     # "Overall_heat_transfer_coeff"       : self.Overall_heat_transfer_coeff(),
        #     # "Heat_duty"       : self.Heat_duty(),
        #     # "Reserve"       : self.Reserve(),
        #     # "Pressure_drop_port"  : self.Pressure_drop_port(),
        #     # "Friction_factor"  : self.Friction_factor(),
        #     # "Pressure_drop_channel"  : self.Pressure_drop_channel(),
        #     # "Total_Pressure_drop"     : self.Total_Pressure_drop(),
        #     # "Capacity"       : self.Capacity(),
        #     # "NTU_Method"        : self.NTU_Method(),
        #     # "get_output"        : self.get_output()
        # }
        # return 


    # def MassFlow(self):
    #     m_prim = self.Heat_load/(self.Heat_cap_water_prim*(self.Inlet_Temp_prim-self.Outlet_Temp_prim))
    #     m_sec = self.Heat_load/(self.Heat_cap_water_prim*(self.Inlet_Temp_sec-self.Outlet_Temp_sec))
    #     return m_prim, m_sec
    
    # def dTLm(self):
    #     dT1 = self.Inlet_Temp_prim-self.Outlet_Temp_sec #[K]
    #     dT2 = self.Outlet_Temp_prim - self.Inlet_Temp_sec #[K]
    #     self.dTlm = (dT1-dT2)/math.log(dT1/dT2)
    #     return dT1, dT2, self.dTlm
    
    # def Eff_No_Plates(self):
    #     Ne = self.Tot_No_Plates - 2
    #     return Ne
    
    # def Pitch(self):
    #     p = (self.Compressed_plate_length/self.Tot_No_Plates)
    #     return p    
        
    # def Mean_channel_spacing(self, p):
    #     b = self.p - self.Plate_Thickness
    #     return b  

    # def Channel_flow_area(self, b):
    #     A_ch = b * self.Horozontal_plate_width
    #     return A_ch   

    # def Single_plate_heat_transfer_area(self, Ne):
    #     A_1 = self.A_tot_eff / Ne
    #     return A_1
    
    # def Projected_plate_area(self):
    #     Lp = self.Vertical_port_length - self.Port_diameter
    #     A_1p = Lp * self.Horozontal_plate_width
    #     return A_1p, Lp
    
    # def Enlargment_factor(self, A_1, A_1p):
    #     enl_factor = A_1 / A_1p
    #     return enl_factor

    # def Hydraulic_diameter_flow_channel(self, b, enl_factor):
    #     Dh = (2 * b) / enl_factor
    #     return Dh
    
    # def Number_channels_per_pass(self):
    #     N_cp = (self.Tot_No_Plates - 1) / (2* self.No_of_passes)
    #     return N_cp
       
    # def Channel_velocity(self):
    #     u_ch_h = self.m_prim / (self.Den_water_prim*self.A_ch*self.N_cp)
    #     u_ch_c = self.m_sec / (self.Den_water_sec*self.A_ch*self.N_cp)
    #     return u_ch_h, u_ch_c

    # def Mass_velocity_per_channel(self, m_prim, N_cp, A_ch, m_sec):
    #     G_ch_h = m_prim / (N_cp * A_ch)
    #     G_ch_c = m_sec / (N_cp * A_ch)
    #     return G_ch_h, G_ch_c

    # def Mass_velocity_port(self, m_prim, m_sec):
    #     A_port = ((math.pi/4) * self.Port_diameter**2)
    #     G_port_h = m_prim / A_port
    #     G_port_c = m_sec / A_port
    #     return G_port_h, G_port_c, A_port

    # def Port_velocity(self, m_prim, m_sec, A_port):
    #     u_p_h = m_prim / (self.Den_water_prim*A_port)
    #     u_p_c = m_sec / (self.Den_water_sec*A_port)
    #     return u_p_h, u_p_c

    # def Prandl(self):
    #     Pr_h = (self.Heat_cap_water_prim * self.Kin_Vis_water_prim * self.Den_water_prim) / self.thermal_cond_water_prim
    #     Pr_c = (self.Heat_cap_water_sec * self.Kin_Vis_water_sec * self.Den_water_sec) / self.thermal_cond_water_sec
    #     return Pr_h, Pr_c

    # def Reynolds(self, u_ch_h, Dh, u_ch_c):
    #     Re_h = (u_ch_h * Dh) / self.Kin_Vis_water_prim
    #     Re_c = (u_ch_c * Dh) / self.Kin_Vis_water_sec
    #     return Re_h, Re_c

    # def Nusselt(self, Re_h, Re_c, Pr_h, Pr_c):
    #     Nu_h = 0.36* (Re_h**(2/3)) * (Pr_h**(1/3))
    #     Nu_c = 0.36* (Re_c**(2/3)) * (Pr_c**(1/3))
    #     return Nu_h, Nu_c

    # def Heat_transfer_coeff(self, Nu_h, Nu_c, Dh):
    #     alpha_h=(Nu_h*self.thermal_cond_water_prim)/ Dh
    #     alpha_c=(Nu_c*self.thermal_cond_water_sec)/ Dh
    #     return alpha_h, alpha_c

    # def Overall_heat_transfer_coeff(self, alpha_h, alpha_c, Rdi, Rdo):
    #     Uo = 1.0/((1.0/alpha_h)+(1.0/alpha_c)+(self.Plate_Thickness/self.Plate_thermal_cond_SS))
    #     Uf= 1.0/((1.0/alpha_h)+(1.0/alpha_c)+(self.Plate_Thickness/self.Plate_thermal_cond_SS)+ Rdi + Rdo)
    #     return Uo, Uf

    # def Heat_duty(self, Uo, Uf, dTlm):
    #     Qc = Uo * self.A_tot_eff * dTlm
    #     Qf = Uf * self.A_tot_eff * dTlm 
    #     return Qc, Qf

    # def Reserve(self, Qc, Qf):
    #     Res_Q = ((Qc/self.Heat_load)-1)*100 
    #     return Res_Q  

    # def Pressure_drop_port(self, G_port_h, G_port_c):
    #     P_p_h = (1.4* self.No_of_passes * self.Den_water_prim * G_port_h**2)  / 2
    #     P_p_c = (1.4* self.No_of_passes * self.Den_water_sec * G_port_c**2)  / 2 
    #     return P_p_h, P_p_c

    # def Friction_factor(self, Re_h, Re_c): #Kumar Korelation
    #     if Re_h < 100:
    #         f_h = 19.400 / (Re_h**0.589)
    #         f_c = 19.400 / (Re_c**0.589)
    #     else:
    #         f_h = 2.99 / (Re_h**0.183)
    #         f_c = 2.99 / (Re_c**0.183)
    #     return f_h, f_c            
        
    # def Pressure_drop_channel(self, f_h, f_c, Dh, G_ch_h, G_ch_c):
    #     P_c_h = 4 * f_h * ((self.Vertical_port_length * self.No_of_passes) / Dh) * (G_ch_h **2 / (2 * self.Den_water_prim)) * (self.Kin_Vis_water_wall_prim / self.Kin_Vis_water_prim)**-0.17 
    #     P_c_c = 4 * f_c * ((self.Vertical_port_length * self.No_of_passes) / Dh) * (G_ch_c **2 / (2 * self.Den_water_sek)) * (self.Kin_Vis_water_wall_sek / self.Kin_Vis_water_sek)**-0.17 
    #     return P_c_h, P_c_c

    # def Total_Pressure_drop(self, P_p_h, P_c_h, P_p_c, P_c_c):
    #     P_t_h = P_p_h + P_c_h
    #     P_t_c = P_p_c + P_c_c
    #     return P_t_h, P_t_c
    
    # def Capacity(self, m_prim, m_sec, P_p_c, P_c_c):
    #     C_h = self.Heat_cap_water_prim*m_prim
    #     C_c = self.Heat_cap_water_sec*m_sec
    #     C_min = min(C_h, C_c)
    #     C_max = max(C_h, C_c)
    #     return C_h, C_c, C_min, C_max

    # def NTU_Method(self, Uo, m_sec, Cmin, Cmax, C_h, C_c, fc):
    #     NTU = (Uo*self.A_tot_eff) / (Cmin)
    #     R = Cmin / Cmax
    #     if self.flow_condition =='counterflow':
    #         P_cc = (1-(math.exp(-NTU*(1-R)))) / (1-(R*(math.exp(-NTU*(1-R))))) #for R:=1 counterflow
    #         Q_new = P_cc * Cmin * (self.Inlet_Temp_prim - self.Inlet_Temp_sec) 
    #     else:
    #         P_co = (1-math.exp(-NTU(1+R))) / (1+R) #for R:=1 co-current
    #         Q_new = P_co * Cmin * (self.Inlet_Temp_prim - self.Inlet_Temp_sec) 
        
    #     T2_new = self.Inlet_Temp_prim - (Q_new / C_h) 
    #     t2_new = self.Inlet_Temp_sec + (Q_new / C_c)
    #     return NTU, R, P_cc, P_co, Q_new, T2_new, t2_new
    
    # def get_output(self):
    #     output = {
    #         "Fouling"  : self.Fouling(),
    #         "MassFlow"  : self.MassFlow(),
    #         "dTLm"     : self.dTLm(),
    #         "Eff_No_Plates"       : self.Eff_No_Plates(),
    #         "Pitch"        : self.Pitch(),
    #         "Mean_channel_spacing"       : self.Mean_channel_spacing(),
    #         "Channel_flow_area"       : self.Channel_flow_area(),     
    #         "Single_plate_heat_transfer_area"  : self.Single_plate_heat_transfer_area(),
    #         "Projected_plate_area"  : self.Projected_plate_area(),
    #         "Enlargment_factor"     : self.Enlargment_factor(),
    #         "Hydraulic_diameter_flow_channel"       : self.Hydraulic_diameter_flow_channel(),
    #         "Number_channels_per_pass"        : self.Number_channels_per_pass(),
    #         "Channel_velocity"       : self.Channel_velocity(),
    #         "Mass_velocity_per_channel"       : self.Mass_velocity_per_channel,
    #         "Mass_velocity_port"  : self.Mass_velocity_port(),
    #         "Port_velocity"  : self.Port_velocity(),
    #         "Prandl"     : self.Prandl(),
    #         "Reynolds"       : self.Reynolds(),
    #         "Nusselt"        : self.Nusselt(),
    #         "Heat_transfer_coeff"       : self.Heat_transfer_coeff(),
    #         "Overall_heat_transfer_coeff"       : self.Overall_heat_transfer_coeff(),
    #         "Heat_duty"       : self.Heat_duty(),
    #         "Reserve"       : self.Reserve(),
    #         "Pressure_drop_port"  : self.Pressure_drop_port(),
    #         "Friction_factor"  : self.Friction_factor(),
    #         "Pressure_drop_channel"  : self.Pressure_drop_channel(),
    #         "Total_Pressure_drop"     : self.Total_Pressure_drop(),
    #         "Capacity"       : self.Capacity(),
    #         "NTU_Method"        : self.NTU_Method(),
    #         "get_output"        : self.get_output()
    #     }
    #     return output        
        
        
        
        
        
        
        
        
        
        
               
        
        
        
        





        
     
        
     
        
     
        
     
        
    

       
