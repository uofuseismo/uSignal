#!/usr/bin/env python3
from scipy.signal import sosfilt
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
    sos_coeffs = [ [0.000401587491686, 0.000803175141692, 0.000401587491549, 1.000000000000000, -1.488513049541281,  0.562472929601870],
                   [1.000000000000000, -2.000000394412897,  0.999999999730209,  1.000000000000000, -1.704970593447777,  0.792206889942566],
                   [1.000000000000000,  1.999999605765104,  1.000000000341065,  1.000000000000000, -1.994269533089365,  0.994278822534674],
                   [1.000000000000000, -1.999999605588274,  1.000000000269794,  1.000000000000000, -1.997472946622339,  0.997483252685326] ]
    x = load_signal('gse2.txt')
    start = time.time_ns()
    y = sosfilt(sos_coeffs, x)
    end = time.time_ns()
    print("Processing time {}".format( (end - start)*1.e-9 ))
    ofl = open('sosReference.gse2.txt', 'w')
    for i in range(len(y)):
        ofl.write("{}\n".format(y[i]))
    ofl.close()
     
