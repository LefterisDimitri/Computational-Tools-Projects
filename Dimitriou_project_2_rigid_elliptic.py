import scipy.special as sp
import numpy as np 
import math as m
import matplotlib.pyplot as plt 

# ---------- STHN ASKHSH AUTH THA LUSOUME TO PROBLHMA TOU STEREOU SWMATOS ME XRHSH TWN ELLEIPTIKWN SUNARTHSEWN JACOBI ---------- #


# orismos metablhtwn kai fusikwn statherwn

# ----- gia to xroniko bhma ----- #
tmax = 100
tmatrix = np.arange(1,tmax+1,1)  # to xroniko vhma to orizw monada, gia mikrotero vhma pairnw kalutera apotelesmata
#print(tmatrix)
#print(len(tmatrix))
#print(tmatrix[99])
# ------------------------------- #

omega_1_t = np.zeros(len(tmatrix)) # arxikopoihsh pinaka xronoeksartomenhs sunistwsas ---> omega_1_t
omega_2_t = np.zeros(len(tmatrix)) # arxikopoihsh pinaka xronoeksartomenhs sunistwsas ---> omega_2_t
omega_3_t = np.zeros(len(tmatrix)) # arxikopoihsh pinaka xronoeksartomenhs sunistwsas ---> omega_3_t

#print(omega_1_t)
#print(omega_1_t[0])
#print(omega_1_t[99])

E_t = np.zeros(len(tmatrix))       # arxikopoihsh tou pinaka gia ton upologismo ths energeias se kathe xroniko bhma
sfalma = np.zeros(len(tmatrix))    # arxikopoihsh tou pinaka sfalmatos
sfalma_t_0 = float("NaN")          # arxikopoihsh tou sfalmatos gia t = 0 ws "kenh timh"

I_1 = 0.8 # roph adraneias ---> I1
I_2 = 0.9 # roph adraneias ---> I2
I_3 = 1.0 # roph adraneias ---> I3

file = open('data_project_2_a.txt','w')
omega_1_t_0 = 1.0 # arxikopoihsh gwniakhs taxuthtas ws pros ton aksona 1 gia thn xronikh stigmh ---> t = 0
omega_2_t_0 = 0.0 # arxikopoihsh gwniakhs taxuthtas ws pros ton aksona 2 gia thn xronikh stigmh ---> t = 0
omega_3_t_0 = 2.0 # arxikopoihsh gwniakhs taxuthtas ws pros ton aksona 3 gia thn xronikh stigmh ---> t = 0

E = 0.5*(I_1*omega_1_t_0**2 + I_2*omega_2_t_0**2 + I_3*omega_3_t_0**2)             # theoritika statherh energeia - arxikopoihsh gia t = 0
M = np.sqrt((I_1*omega_1_t_0)**2 + (I_2*omega_2_t_0)**2 + (I_3*omega_3_t_0)**2)    # theoritika statherh stroformh - arxikopoihsh gia t = 0
file.write('%d %f %f %f %f\n'%(0, omega_1_t_0, omega_2_t_0, omega_3_t_0, sfalma_t_0)) 

for i in range(0, len(tmatrix)): 

    # ---------- Jacobi elliptic functions ---------- #
    tau = tmatrix[i]*np.sqrt((I_3-I_2)*(M**2-2.0*E*I_1)/(I_1*I_2*I_3))
    k_sq = (I_2-I_1)*(2.0*E*I_3-M**2)/((I_3-I_2)*(M**2-2.0*E*I_1))
    [s, c, d, phi] = sp.ellipj(tau, k_sq)
    omega_1_t[i] = np.sqrt((2*E*I_3 - M**2)/(I_1*(I_3 - I_1)))*c
    omega_2_t[i] = np.sqrt((2*E*I_3 - M**2)/(I_2*(I_3 - I_2)))*s
    omega_3_t[i] = np.sqrt((M**2 - 2*E*I_1)/(I_3*(I_3 - I_1)))*d
    # ----------------------------------------------- #
    
    E_t[i] = 0.5*(I_1*omega_1_t[i]**2 + I_2*omega_2_t[i]**2 + I_3*omega_3_t[i]**2)    # ypologismos ths energeias
    
    # ---------- gia to sfalma ---------- #
    d = E_t[i]-E
    if d == 0:
        sfalma[i] = float("NaN")
       
    else:
        sfalma[i] = m.log10(abs((E_t[i]-E)/E))
    # ----------------------------------- #    
    
    file.write('%d %f %f %f %f\n'%(i, omega_1_t[i], omega_2_t[i], omega_3_t[i], sfalma[i])) # arxeio me 101 times ( t = 0,...,100 )
   
file.close()

data = np.loadtxt('data_project_2_a.txt')
x = data[:,0]
y1 = data[:,1]
plt.figure(1)
plt.plot(x,y1)
plt.xlabel("t")
plt.ylabel("ω1")
plt.show()

y2 = data[:,2]
plt.figure(2)
plt.plot(x,y2)
plt.xlabel("t")
plt.ylabel("ω2")
plt.show()

y3 = data[:,3]
plt.figure(3)
plt.plot(x,y3)
plt.xlabel("t")
plt.ylabel("ω3")
plt.show()
    
y4 = data[:,4]
plt.figure(4)
plt.plot(x,y4)
plt.xlabel("t")
plt.ylabel("error")
plt.show()
    
#plt.plot(tmatrix,omega_1_t)
#plt.show()

#plt.plot(tmatrix,omega_2_t)
#plt.show()

#plt.plot(tmatrix,omega_3_t)
#plt.show()

#plt.plot(tmatrix,sfalma)
#plt.show()
