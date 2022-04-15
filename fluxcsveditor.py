#Fenfang Wu, GPS Caltech, 2020
#!env python3
import numpy as np
from numpy import *
from tqdm import tqdm
from scipy import optimize
import pylab as py
from scipy import integrate
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import csv
from numpy import zeros,array,log,concatenate
from copy import deepcopy

def fluxconvert(kf,kr,networkbuilder):
    new_rows = []
    count = 0
    for i in networkbuilder:
        if i == '.':
            value = count
        count += 1
    
    filename = networkbuilder[0:value]
    networkbuilt = filename + '.csv'
    with open(networkbuilt,newline='') as f:
        reader = csv.reader(f) # pass the file to our csv reader
        counter = 0
        for row in reader:     # iterate over the rows in the file
            new_row = row
            if row[0] != 'Reaction':

                    new_row[1] = kf[counter]
                    new_row[2] = kr[counter]
                    counter += 1
            
            new_rows.append(new_row) # add the modified rows

    newfile = filename + '_FLUXINVERTED.csv'
    with open(newfile, 'w') as f:
        # Overwrite the old file with the modified rows
        writer = csv.writer(f)
        writer.writerows(new_rows)
