#!/usr/bin/python
# python code using mpmath to create a test data sets
from mpmath import *
mp.dps = 30
mp.pretty = True

import random
random.seed(1)

N = 1000
mu = 0.0
sigma = 2.0
fid = open('li_testdata.dat', 'w')
fid.write('# test data from li_testdata.py\n')
fid.write('# using mpmath library\n')
fid.write('# s(real, imag), z(real, imag), Li(real, imag)\n')
for i in range(1,N):
    s = random.gauss(mu, sigma) + j*random.gauss(mu, sigma)
    z = random.gauss(mu, sigma) + j*random.gauss(mu, sigma)
    L = polylog(s, z)
    fid.write('%.20f,%.20f,%.20f,%.20f,%.20f,%.20f\n' % (s.real, s.imag, z.real, z.imag, L.real, L.imag))


fid.close()


