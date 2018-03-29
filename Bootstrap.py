# -*- coding: utf-8 -*-
"""
Created on Sat Mar 10 20:05:15 2018

@author: alexa

Bootstrap function takes in csv datafile with no header, and B which is the number of replicates wanted.


"""
import numpy as np
from scipy import genfromtxt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



def Bootstrap(data, B):
    N, D = np.shape(data)
    
    #Set up training set matrix for bootstrapping
    BTRAIN = np.zeros((B,D))
    
    #Create Matrix for determining the center of the data
    CNTER = np.zeros((D,1))
    
    #Create Replicate matrix
    BSAMP = np.zeros((N,D))
    
    #Determine picks of the dataset for each replicate
    PICKS = np.random.randint(0,N,(B,N))
    
    #For each replicate
    for i in range(0,B):
        #Take in the data from our picks from the current replicate to end of the picks on
        BSAMP = data[PICKS[i],:]
        #Sum these numbers on the column axis, and then average them by the number of spectra given.
        BTRAIN[i] = np.sum(BSAMP,axis=0)/N
    
    #Determine the Center by Averaging the Columns of the new Replicate Data.
    for i in range(len(CNTER)):
        CNTER = np.sum(BTRAIN,axis=0)/B
    #Return the     
    return BTRAIN, CNTER, BSAMP



