#!/usr/bin/env python3
from scipy.signal import lfilter
from numpy import zeros
import time

def load_signal(file_name = 'gse2.txt'):
    ifl = open(file_name, 'r')
    data = ifl.read()
    data = data.split('\n')
    signal = zeros(len(data) - 1)
    for i in range(len(signal)):
        signal[i] = float(data[i])
    return signal
    
if __name__ == "__main__":
    b = [0.000401587491686, 0.0, -0.001606349966746, 0.0,
         0.002409524950119, 0.0, -0.001606349966746, 0.0, 
         0.000401587491686];
    a = [1.000000000000000, -7.185226122700763, 22.615376628798678, -40.733465892344896,
         45.926605646620146,-33.196326377161412, 15.023103545324197, -3.891997997268024, 
         0.441930568732716]
    x = load_signal('gse2.txt')
    start = time.time()
    y = lfilter(b, a, x)
    end = time.time()
    ofl = open('transposeDF2.gse2.txt', 'w')
    for i in range(len(y)):
        ofl.write("{}\n".format(y[i]))
    ofl.close()
     
