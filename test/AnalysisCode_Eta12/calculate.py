#!/bin/env python

import os, sys, string, math

k = string.atoi(sys.argv[1])
n = 25.0
p = 0.4

def binomial(n,k):
   toReturn = math.factorial(n)/math.factorial(k)
   return toReturn/math.factorial(n - k)

def calculate(m):
   return binomial(25, m)*(p**m)*((1 - p)**(n-m))

toReturn = 0.0
for m in range(0,k):
   toReturn += calculate(m)

print 1.0-toReturn
