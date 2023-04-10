# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 12:32:53 2023

@author: kieran
"""
import numpy as np
from matplotlib import pyplot as plt

def phiupdate(phidot, phi, DelT):
    phinew = DelT*phidot + phi
    return(phinew)

def phidotupdate(H, phidot, phi, DelT):
    Vphi = -2*phi**-3
    phidotnew = phidot -(3*H*phidot + (20)*Vphi)*DelT
    return(phidotnew)

def A_update(H, a, DelT):
    a_new = DelT*a*H + a
    return(a_new)

def Friedmann(H0, omr, omm, oml, a):
    Hsq = H0**2*(omr*(a**-4) + omm*(a**-3) + oml)
    return(Hsq**0.5)

def density(phidot, phi):
    Vphi = phi**-2
    density = 0.5*(phidot**2) + Vphi*(20)
    return(density)

def VVphi(phi):
    VVphis = phi**-2
    return(VVphis)

def main():
    
    #Friedmann Parameters
    H0 = 0.07 #For testing the integrator, set H0 as 1
    omr = 0
    omm = 0.3149
    oml = 0.6851 #Set oml as 1
    a = 0.001 #Set a as 1
    DelT = 0.01
    a_new = a
    
    #general parameters
    t = 0
    i = 0
    
    #Klein Gordon Parameters
    #Initial condition of 4.3*10^-4
    phi0 = 1
    phidot0 = 0
    VV0 = VVphi(phi0)
    phiden = density(phidot0, phi0)
    
    #sets oml properly
    oml = (8/3)*np.pi*phiden
    print(str(oml) + "oml")
    
    #sets omm properly
    matter_par = 0.3149*(a**-3)
    matter_density = matter_par/((8/3)*np.pi)

    
    
###############################################################################
#Setting up lists
###############################################################################
    Hubble_List = [H0]
    Time_List = [t]
    a_list = [a]
    phi_list = [phi0]
    phidot_list = [phidot0]
    Potential_List = [VV0]
    phiden_list = [phiden]
    matter_density_list = [matter_density]
    
###############################################################################
#INITIAL UPDATE
###############################################################################
    #updates friedmann parameters
    H = Friedmann(H0, omr, omm, oml, a)
    print(str(H) + " " + "H")
    a = A_update(H, a, DelT)
    print(str(a) + " " + "a")
    
    #updates klein gordon parameters
    phi = phiupdate(phidot0, phi0, DelT)
    Potential_List.append(VVphi(phi))
    print(str(phi) + " " + "phi")
    phidot = phidotupdate(H, phidot0, phi, DelT)
    print(str(phidot) + " " + "phidot")
    t += DelT
    
    #updates dark energy densities
    phiden = density(phidot, phi)
    print(str(phiden) + "phiden")
    oml = (8/3)*np.pi*phiden
    print(str(oml) + "oml")
    
    #updating matter densities
    matter_par = 0.3149*(a**-3)
    matter_density = matter_par/((8/3)*np.pi)
    
    #updates relevant lists
    Hubble_List.append(H)
    a_list.append(a)
    phi_list.append(phi0)
    phidot_list.append(phidot0)
    Time_List.append(t)
    phiden_list.append(phiden)
    matter_density_list.append(matter_density)

    
###############################################################################
#Integrator
###############################################################################
    while i <= 2500: #runs for 20 timesteps for now
        
        #friedmann updates
        H = Friedmann(H0, omr, omm, oml, a) #calculates new H
        print(str(H) + " " + "H") 
        a = A_update(H, a, DelT) #updates A
        print(str(a) + " " + "a")
        
        #update klein gordon
        phi = phiupdate(phidot, phi, DelT)
        Potential_List.append(VVphi(phi))
        phidot = phidotupdate(H, phidot, phi, DelT)
        print(str(phidot) + " " +  "phidot")
        print(str(phi) + " " + "phi")
        
        #updates time
        t += DelT 
        
        #dark energy density update
        phiden = density(phidot, phi)
        oml = (8/3)*np.pi*phiden
        print(str(oml) + " " + "oml")
        
        #updating matter densities
        matter_par = 0.3149*(a**-3)
        matter_density = matter_par/((8/3)*np.pi)
        
        #appends_lists
        Hubble_List.append(H)
        Time_List.append(t)
        a_list.append(a)
        phi_list.append(phi)
        phidot_list.append(phidot)
        phiden_list.append(phiden)
        matter_density_list.append(matter_density)
        
        #Moves onto next iteration
        i += 1
        
        #Plotting potential against time
        plt.plot(Time_List, Potential_List)
        plt.show()
        plt.cla
        
    matter_array = np.array(matter_density_list)
    phiden_array = np.array(phiden_list)
    a_array = np.array(a_list)
    
###############################################################################
#Plots
###############################################################################

    plt.plot(Time_List, a_list)
    plt.title("Plot of scale factor varying with time")
    plt.xlabel("Time(GY)")
    plt.ylabel("Scale Factor")
    plt.show()
    
    
    plt.plot(Time_List, matter_density_list)
    plt.title("Plot of matter density varying with time")
    plt.xlabel("Time(GY)")
    plt.ylabel("matter density")
    plt.show()
    
    plt.plot(Time_List, phiden_list)
    plt.title("Plot of DE density varying with time")
    plt.xlabel("Time(GY)")
    plt.ylabel("DE density")
    plt.show()
    
    plt.plot(np.log(matter_array), np.log(phiden_array))
    plt.xlabel("matter density")
    plt.ylabel("DE density")
    plt.show()
    
    plt.plot(np.log(a_array), np.log(matter_array), label="Matter Density")
    plt.plot(np.log(a_array), np.log(phiden_array), label="Dark Energy Density")
    plt.title("Plot of Densities against Scale Factor")
    plt.xlabel("Density")
    plt.ylabel("Scale Factor")
    plt.legend(loc="upper right")
    plt.show()
    
    print(a_array)
    print(np.log(a_array))
    print(phiden_array)
    
    
main()







