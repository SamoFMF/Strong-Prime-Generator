# Strong Prime Generator

from random import randrange
from math import log, ceil, sqrt

class SPG:

    @staticmethod
    def MR(n, K=1, rand=randrange):
        '''Miller-Rabin Algorithm.'''

        # Hardcoded list of <=10-bit primes
        if n > 1024:
            ps = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021]
            for p in ps:
                if n%p == 0:
                    return False
        m = n-1
        k = 0
        while m%2 == 0:
            m //= 2
            k += 1
        for _ in range(K):
            a = rand(1,n)
            b = pow(a,m,n)
            if b == 1:
                continue
            for _ in range(k):
                if b == n-1:
                    break
                else:
                    b = pow(b,2,n)
            else:
                return False
        return True

    @staticmethod
    def generateStrongPrime(n=512, pon1=20, pon2=20, rand=randrange):
        '''Generates a n-bit strong prime.'''
        
        log2 = lambda x: ceil(log(x,2)) # base 2 log

        n1 = (n - log2(n)) // 2 - 4 # n1+n1+1 < n - log2n - 6
        n2 = n1 - log2(n1) - 7 # 6 to be sure we find r of desired bit length

        # Find n1-bit prime s
        M = 2**(n1-2)
        M2 = 2*M
        a = rand(M2, M+M2) # interval [2**(n-1), 3*2**(n-2)) | 3*2**(n-2) = 2**(n-1) + 2**(n-2)
        s = a+1 if a%2 == 0 else a+2 # s is a "large" prime => must be odd
        M2 *= 2 # The upper bound: 2**n1
        while s < M2:
            if SPG.MR(s,pon1):
                break
            s += 2 # Skip even values
        else:
            # We made sure that this should never happen - but you never know?
            raise RuntimeError("Generating prime s exceeded bit size.")
        
        # Find n2-bit prime t
        M = 2**(n2-2)
        M2 = 2*M
        b = rand(M2, M+M2) # interval: same as above for s - sub n1 for n2
        t = b+1 if b%2 == 0 else b+2 # t is a "large" prime => must be odd
        M2 *= 2 # The upper bound: 2**n2
        while t < M2:
            if SPG.MR(t,pon1):
                break
            t += 2 # Skip even values
        else:
            # We made sure that this should never happen - but you never know?
            raise RuntimeError("Generating prime t exceeded bit size.")
        
        # Construct r
        # Idea described in my paper
        r = 1
        rlow = r
        rhigh = r
        blow = 2**(n1-1)
        bhigh = 2*blow
        m = 0
        step = n1-n2-1
        t2 = t * 2**(step+1)
        while step > 0:
            if rlow+t2 < blow:
                rlow += t2
                m += 2**step
            if rhigh+t2 < bhigh:
                rhigh += t2
            step -= 1
            t2 //= 2
        while rlow < blow:
            rlow += t2
            m += 1
        while rhigh < bhigh:
            rhigh += t2
        r1 = m*t2 + 1
        while not SPG.MR(r1, pon1):
            r1 += t2
            if r1 > bhigh:
                raise RuntimeError("Generating prime r exceeded bit size.")
        r = r1
        r1 = 0
        
        # Construct n-bit prime p from r and s
        rs = r*s
        urs = (pow(s,r-1,rs) - pow(r,s-1,rs)) % rs
        if urs%2 == 0:
            p0 = urs + rs
        else:
            p0 = urs
        rs *= 2
        p = p0
        M2 = 2**(n-1)

        plow = p
        phigh = p
        if n > 65:
            # For n around 1.5k we start running into problems with large floats
            # That's why we seperate sqrt(2) * 2**(n-1) into two parts
            # sqrt(2) * 2**64 * 2**(n-1-64) and use "ceil" on the first
            # two elements - this leaves us with 2 integers and we avoid
            # having to deal with large floating point values
            blow = ceil(sqrt(2) * 2**64) # the error is 0% on this
            blow *= 2**(n-65) # lower boundary
        else:
            blow = sqrt(2) * 2**(n-1) # lower boundary
        bhigh = 2**n # upper boundary
        m = 0 # lower boundary multiplier
        step = n-2*n1
        rs *= 2**step
        while step > 0:
            # Lower boundary
            if plow+rs < blow:
                plow += rs
                m += 2**step
            # Upper boundary
            if phigh+rs < bhigh:
                phigh += rs
            step -= 1
            rs //= 2
        while plow < blow:
            plow += rs
            m += 1
        while phigh < bhigh:
            phigh += rs
        
        p1 = p + m*rs
        while not SPG.MR(p1,pon2):
            p1 += rs
            if p1 > bhigh:
                raise RuntimeError("Generating prime p exceeded bit size.")
        p = p1
        p1 = 0

        return p,r,s,t