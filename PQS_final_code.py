import numpy as np
import math

#    *************************  PQS | Assignment 1 | Souvik | 231110405   *************************
n = 512
q = 12289  #(Z 12289)
r = 10968  #gamma

br=[] #bitreverse array

p1 = np.random.randint(0, q, n)
p2 = np.random.randint(0, q, n)

#    *************************   Negative-wrapped Convolution (Verification) *************************
f = np.zeros(n + 1)
f[0] = 1
f[n] = 1
ans = np.remainder(np.polydiv(np.polymul(p1[::-1],p2[::-1]), f)[1], q).astype(int)[::-1]

print()
print("   ********************     Algorithm Starts     ********************   ")
print("")
print("Initial 'P1': ",p1)
print("")
print("and Initial 'P2': ",p2)

gamma=[]                                                #Gamma matrix
for i in range(n):
    gamma.append((r**i)%q)
print("\n")
print("The Gamma array: ",gamma)

def BiReA(po):                                          #bit reverse sequence 
    df=[]
    dfd=[]
    for i in po:
        df.append(i)
    for i in range(n):                                         
        s=str(bin(i)[2:])
        g=int(math.log(n,2))-len(s)
        for j in range(g):
            s="0"+s
        s=s[::-1]
        s=int(s,2)
        dfd.append(df[s])    
    return(dfd) 

gamma=BiReA(gamma)
print("")
print("The Gamma array after Bit Reverse is: ",gamma)

def NTT(pq):                                            #DIF-NTT Transformation code
    t=n
    m=1
    while(m<n):
        t=t//2
        for i in range(m):
            j1=2*i*t
            j2=j1+t-1
            s=gamma[m+i]
            for j in range(j1,j2+1):
                u=pq[j]
                v=pq[j+t]*s
                pq[j]=(u+v)%q
                pq[j+t]=(u-v)%q
        m=2*m
    return(pq)

x=NTT(p1)
y=NTT(p2)

print("")
print("'P1' after NTT: ",x)
print("")
print("'P2' after NTT: ",y)

#    *************************   Point wise multiplication   *************************

def PWM(p11,p12):                                            #Point wise multiplication
    pw=[]
    for i in range(len(p11)):
        pw.append((p11[i]*p12[i])%q)
    return(pw)

p=PWM(x,y)
print("")
print("The array after 'NTT of P1 & P2', and 'Point Wise Multiplication of NTT' : ",p)

#    *************************   DIF-INTT    *************************

# gamma inverse (r*r_inv)mod q=1
r_inv=1
for i in range(2,q):                                       #finding omega inverse
    if(i*r)%q==1:
        r_inv=i 
        exit

# n inverse (n*n_inv)mod q=1
n_inv=1
for i in range(2,q):                                       #finding omega inverse
    if(i*n)%q==1:
        n_inv=i 
        exit

gammainv=[]  
for i in range(n):                                               #Gamma Inverse matrix
    gammainv.append((r_inv**i)%q)

print("")
print("The Gamma Inverse array: ",gammainv)
gammainv=BiReA(gammainv)
print("")
print("'Gamma Inverse Array' after Bit Reverse: ",gammainv)

def INTT(po):                                       #DIF, I-NTT Transformation code
    t=1
    m=n
    while(m>1):
        j1=0
        h=m//2
        for i in range(h):
            j2=j1+t-1
            s=gammainv[h+i]
            for j in range(j1,j2+1):
                u=po[j]
                v=po[j+t]
                po[j]=(u+v)%q
                po[j+t]=((u-v)*s)%q
            j1=j1+(2*t)
        t=2*t
        m=m//2
    for i in range(n):
        po[i]=(po[i]*n_inv)%q
    return(po)    

final=INTT(p)

print("")
print("'Final array, after 'NTT of P1 & P2', 'Point Wise Multiplication of NTT', and 'Inverse-NTT of the o/p from NWC':",final)
print("")

print("After negative wrapping 'i.e. after Modulo(X^n+1)' (Actual Answer Array for Verification): ",ans)
print("")
print(" **** Verification, whether \"NTT, PWM and INTT\" is same to \"Negative Wrapped Convolution **** ")
print("Both are same: ",np.array_equal(final, ans))


