# RSA Using Strong Prime Generator

from spg import SPG
from random import randrange
from math import gcd

class RSA:

    def __init__(self, n, pon1=20, pon2=10, rand=randrange):
        while True:
            self.e = rand(2**16+1, 2**256)
            if self.e % 2 == 1: break
        self.n = n
        self.pon1 = pon1
        self.pon2 = pon2
        self.p = None
        self.q = None
        self.key = None
        self.phi = None
    
    def getPrimes(self):
        while True:
            p,_,_,_ = SPG.generateStrongPrime()
            if gcd(p-1,self.e) == 1:
                self.p = p
                break
        while True:
            q,_,_,_ = SPG.generateStrongPrime()
            if abs(self.p-q)>2**(self.n//2 - 100) and gcd(q-1,self.e) == 1:
                self.q = q
                break
        self.key = p*q
        self.phi = (p-1)*(q-1)