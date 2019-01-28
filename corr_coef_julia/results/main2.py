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


def periodogram(t,c,name):
    F=[]
    A=[]
    
    window = np.kaiser(c.shape[-1],5)
    c_win  = c * window


    f,a=DO_FFT(c_win, 1)
    F.append(f)
    A.append(a)

    plt.plot(F[0],A[0],'k', markersize=0.1)
    plt.show()

    #try:
    #    plt.savefig("graph/"+name+"_FFT.png", figuresize=(8,6), dpi=320, format="png")
    #    print_ok("Graph saved in: "+"graph/"+name+"_FFT.png")
    #except IOError as IoE:
    #    print_fail("I/O Error! Erro number = {0}; {1}".format(IoE.errno,IoE.strerror))
    #    exit(IoE.errno)






if __name__ == "__main__":
    
    t,c=load_file(sys.argv[1])
    
    try:
        periodogram(t,c,"corr")
        print ""
        print_ok("Stage1 completed - FFT")
        print ""
    except:
        print_fail("Some eorros has been occured in the Stage1")
        exit(1)


    exit(0)
