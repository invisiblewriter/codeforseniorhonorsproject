# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 17:08:00 2023
Final file: general form of the Friedmann equation using finite differencing methods.
Plotting of the Hubble constant(FINISH TODAY AND THEN SLOW ROLL TOMORROW)

@author: kiera
"""
import numpy as np
import math
from matplotlib import pyplot as plt
import scipy as sp


def Second_Order_Finite_Differencing(anext, anow, aprev, DelT):
    fdashdash = (anext -2*anow + aprev)/(DelT**2)
    return(fdashdash)

def A_update(H, a, DelT):
    """
    Parameters
    ----------
    H : Calculated Value of the Hubble constant
    a : initial scale factor
    DelT : change in time

    Returns
    -------
    Updated scale factor, a(t + DelT) for later use in integrator.
    Add in error message for later

    """
    a_new = DelT*a*H + a
    return(a_new)

def Friedmann(H0, omr, omm, oml, a):
    """
    Parameters
    ----------
    H0 : Initial hubble constant value
    omr : radiation density in universe
    omm : matter/dark matter(I think) density in universe
    oml : dark energy density in universe
    a : scale factor
    Returns
    -------
    A value of the Hubble parameter, determined by the integration

    """
    Hsq = H0**2*(omr*(a**-4) + omm*(a**-3) + oml)
    return(Hsq**0.5)

def DEanalytic(a, t):
    y = (1.31)*t**(2/3) 
    return(y)

def Linear_Fit(x, M):
    yfin = M*x 
    return yfin

def main():

###############################################################################
#This program will plot the how the hubble parameter has changed with time,
#according to the Friedmann equation.

#SETTING INITIAL CONDITIONS/VALUES
###############################################################################
    H0 = 1.0 #For testing the integrator, set H0 as 1
    omr = 0
    omm = 1
    oml = 0 #Set oml as 1
    a = 0.001 #Set a as 1
    DelT = 0.01
    t = 0
    i = 0
    H = H0
    a_new = a
    
    """Criticisms, a can never be zero, this will make everything zero always. H0
    can never be zero as everything would also be zero always. 
    If they end up at zero they'll stay that way so look out for this"""
    """DelT can never be zero also, this will ensure nothing changes""" 
###############################################################################
#Creating lists
###############################################################################
    Hubble_List = [H0]
    Time_List = [t]
    a_list = [a] #This list is for comparisons with analytic solutions
    delA_list = [0]
    ANALYTIC_LIST = [DEanalytic(a, t)]
    a_newnew = a
    #printing lists for debugging
    #print(Hubble_List)
    #print(Time_List)
###############################################################################
#INTEGRATOR
###############################################################################
    while i <= 2000: #runs for 20 timesteps for now
        a_new = A_update(H, a_new, DelT) #updates A
        #print(a_new)
        delArat = (a_new - a_newnew)/a_newnew
        H = Friedmann(H0, omr, omm, oml, a_new) #calculates new H
        t += DelT #updates time
        
        #appends_lists
        Hubble_List.append(H)
        Time_List.append(t)
        a_list.append(a_new)
        delA_list.append(delArat)
        ANALYTIC_LIST.append(DEanalytic(a_new, t))
        a_newnew = a_new
        
        #Moves onto next iteration
        i += 1
###############################################################################
#Plotting
###############################################################################
    
    #printing lists for debugging
    #print(Hubble_List)
    #print(a_list)
    #print(delA_list)
    #print(Time_List)
    
    plt.plot(Time_List, Hubble_List)
    plt.title("Plot of Hubble Parameter Value as varying with time")
    plt.xlabel("Time(UNITS?)")
    plt.ylabel("Hubble parameter(UNITS?)")
    plt.show()
    
    plt.plot(Time_List, a_list, color = "black")
    plt.xlabel("Time(Gy)")
    plt.ylabel("a(t)")
    plt.show()
    
    plt.plot(Time_List, delA_list)
    plt.title("Plot of DelA/A against T")
    plt.xlabel("DelA/A")
    plt.ylabel("Time(UNITS)")
    plt.show()
    
    plt.plot(Time_List, ANALYTIC_LIST)
    plt.plot(Time_List, a_list)
    plt.show()
###############################################################################
#Independence hypothesis testing/fitting
###############################################################################
    #print(ANALYTIC_LIST)
    h = 0
    
    """while h < 19500:
        ANALYTIC_LIST.pop(0)
        h+=1
    h = 0"""

    """while h < 19500:
        a_list.pop(0)
        h+=1
    h = 0"""

    #print(a_list)
    
    delAAarray = np.array(ANALYTIC_LIST)
    a_array = np.array(a_list)
    
    delAarray = delAAarray - a_array
    #print(delAarray)
    percentage_error_array = 100*delAarray/ANALYTIC_LIST
    
    #Removes zero-zero element, as this cannot be analysed. You need to mention
    #that you've removed singularities in your report, when going over error
    #analysis, but say WHY you can ignore them. In this case, a difference between
    #analytic and measured of zero would cause an issue in the error finding functions
    #We can remove this zero as a zero difference does not reflect an infinite error, it
    #reflects a zero error
    """
    while h < 19500:
        Time_List.pop(0)
        h+=1
    h = 0"""
    
    #delAAarray.pop(0)
    
    #print(percentage_error_array)
    
    plt.plot(Time_List, percentage_error_array, color="black")
    plt.ylabel("Percentage Error in a(t)")
    plt.xlabel("Time(Gy)")
    plt.show()
    
    plt.plot(np.log(Time_List), np.log(percentage_error_array), color="black")
    plt.xlabel("ln(T/T0)")
    plt.ylabel("ln(% error)")
    plt.title("Convergence")
    plt.show()
    
    plt.plot(Time_List, delAarray)
    plt.show()
    
###############################################################################
#Code For Determining points of inflection, where new relationships begin to dominate
###############################################################################
    
    K = len(a_list) - 1
    inflection_list = []
    i = 1
    while i < K:
        fne_inflection = Second_Order_Finite_Differencing(a_list[i+1], a_list[i], a_list[i-1], DelT)
        print(fne_inflection)
        inflection_list.append(fne_inflection)
        i += 1
    print(inflection_list)
    #Time_List.pop(0)
    #Time_List.pop(-1)
    
    #Determining inflection points: Write the code. Look at the final value for what
    #actually goes on, is it a true inflection point or an artefact of our data somehow?
    
###############################################################################
#Curve Fitting relationship testing
###############################################################################
    #param, param_cov = sp.optimize.curve_fit(Linear_Fit, Time_List, percentage_error_array)
    #print(param)
    
    plt.plot(Time_List, Hubble_List)
    plt.xlabel("Time(Gy)")
    plt.ylabel("H(1/Gy)")
    plt.show()
    
    #plt.plot(np.log(Time_List), np.log(Hubble_List))
    #plt.xlabel("lg(Time/T0)")
    #plt.ylabel("lg(H/H0)")
    
    Time_List.pop(0)
    Time_List.pop(0)
    Time_List.pop(0)
    Time_List.pop(0)
    Time_List.pop(0)
    Time_List.pop(0)
    Time_List.pop(0)
    Time_List.pop(0)
    
    Hubble_List.pop(0)
    Hubble_List.pop(0)
    Hubble_List.pop(0)
    Hubble_List.pop(0)
    Hubble_List.pop(0)
    Hubble_List.pop(0)
    Hubble_List.pop(0)
    Hubble_List.pop(0)
    
    plt.plot(Time_List, Hubble_List)
    plt.xlabel("Time(Gy)")
    plt.ylabel("H(1/Gy)")
    plt.show()
    
    #plt.plot(np.log(Time_List), np.log(Hubble_List))
    #plt.xlabel("lg(Time/T0)")
    #plt.ylabel("lg(H/H0)")
    #plt.show()
    
    print(Hubble_List[-1])
    
    

main()
    
    
