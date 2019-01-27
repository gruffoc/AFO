#!/usr/bin/python

#default module
import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack as fftp
import scipy.optimize as opt
from scipy import signal
import sys
import os
import string

from bcolors import *
from funct import *
from histogram import *


def periodogram(t,c,ce,name):
    components=[c]
    F=[]
    A=[]

    for i in range(0,len(components)):
        #window=np.kaiser(components[i].shape[-1],50)
        c_coef=components[i]#*window
        f,a=DO_FFT(c_coef,1)
        F.append(f)
        A.append(a)

    plt.plot(t,c_coef)
    plt.show()

    
    plt.figure(1,figsize=(11,7))
    plt.title("Correlation Coefficient $C_{00}^{0t}$ FFT")
    #plt.subplot(121)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Power Spectrum [dB]")
    plt.semilogy(F[0],A[0],'k', markersize=0.1)
    plt.show()


    #try:
    #    plt.savefig("graph/"+name+"_FFT.png", figuresize=(8,6), dpi=320, format="png")
    #    print_ok("Graph saved in: "+"graph/"+name+"_FFT.png")
    #except IOError as IoE:
    #    print_fail("I/O Error! Erro number = {0}; {1}".format(IoE.errno,IoE.strerror))
    #    exit(IoE.errno)






if __name__ == "__main__":
    
    t,c,ce=load_file(sys.argv[1])
    
    c_norm = c

    try:
        periodogram(t,c_norm,ce,"corr_coef")
        print ""
        print_ok("Stage1 completed - FFT")
        print ""
    except:
        print_fail("Some eorros has been occured in the Stage1")
        exit(1)


    exit(0)
