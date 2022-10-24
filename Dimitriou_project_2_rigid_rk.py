import scipy.special as sp
import numpy as np 
import math as m
import matplotlib.pyplot as plt 

# ---------- STHN ASKHSH AUTH THA LUSOUME TO PROBLHMA TOU STEREOU SWMATOS ME XRHSH THS RUNGE-KUTTA 4HS TAKHS ---------- #


# orismos metablhtwn kai fusikwn statherwn ws global metablhtes

I_1 = 0.8 # roph adraneias ---> I1
I_2 = 0.9 # roph adraneias ---> I2
I_3 = 1.0 # roph adraneias ---> I3

tmax = 100
h = 0.1     # to xroniko vhma gia thn RK4 ---> dt == h
N = 1000    # N = tmax/h

sfalma = np.zeros(1000)
E_t = np.zeros(1000)

tmatrix = np.arange(0.1,tmax+0.1,0.1) # xroniko vhma ---> dt = 0.1
#print(tmatrix)
#print(len(tmatrix))

omega_vector = np.zeros(3)
omega_vector_n = np.zeros(3)

#print(len(omega_vector))

file = open('data_project_2_b.txt','w')
omega_vector[0] = 1.0   # arxikopoihsh tou omega_1_t gia thn xronikh stigmh ---> t = 0
omega_vector[1] = 0.0   # arxikopoihsh tou omega_2_t gia thn xronikh stigmh ---> t = 0
omega_vector[2] = 2.0   # arxikopoihsh tou omega_3_t gia thn xronikh stigmh ---> t = 0

E = 0.5*(I_1*omega_vector[0]**2 + I_2*omega_vector[1]**2 + I_3*omega_vector[2]**2)       # energeia gia t = 0 ---> theorhtika statherh

sfalma_t_0 = float("NaN")   # sfalma gia t = 0

file.write('%f %f %f %f %f\n'%(0.0, omega_vector[0], omega_vector[1], omega_vector[2], sfalma_t_0)) 

dt = 0.0               # arxikopoihsh enos xronikou vhmatos gia thn grafikh parastash

def RK(omega_vector, omega_vector_n):   # sunarthsh gia thn RK4
    k1 = np.zeros(3)
    k2 = np.zeros(3)
    k3 = np.zeros(3)
    k4 = np.zeros(3)
    dot_omega_vector = np.zeros(3)
    Y = np.zeros(3)
    
    for i in range(0, 3): 
        Y[i] = omega_vector[i] 
           
    for i in range(0, 3):
        derivatives(Y,dot_omega_vector)
        k1[i] = h*dot_omega_vector[i]
        
    for i in range(0, 3): 
        Y[i] = omega_vector[i] + 0.5*k1[i]
        
    for i in range(0, 3):
        derivatives(Y,dot_omega_vector)
        k2[i] = h*dot_omega_vector[i]
    
    for i in range(0, 3): 
        Y[i] = omega_vector[i] + 0.5*k2[i]
        
    for i in range(0, 3):
        derivatives(Y,dot_omega_vector)
        k3[i] = h*dot_omega_vector[i]
    
    for i in range(0, 3): 
        Y[i] = omega_vector[i] + k3[i]
        
    for i in range(0, 3):
        derivatives(Y,dot_omega_vector)
        k4[i] = h*dot_omega_vector[i]
        
    for i in range(0, 3): 
        omega_vector_n[i] = omega_vector[i] + (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0
    
        
        
def derivatives(Y, dot_omega_vector):              # sunarthsh gia ton upologismo twn paragwgwn
    dot_omega_vector[0] = (I_2-I_3)*Y[1]*Y[2]/I_1
    dot_omega_vector[1] = (I_3-I_1)*Y[2]*Y[0]/I_2
    dot_omega_vector[2] = (I_1-I_2)*Y[0]*Y[1]/I_3
    
for p in range(0, N):
    
    # ----- RK4 ----- # 
    RK(omega_vector, omega_vector_n)
    # --------------- #
    
    # ----- ananewsh arxikwn sunthikwn gia thn epomenh epanalhpsh ----- # 
    for j in range(0, 3):
        omega_vector[j] = omega_vector_n[j] 
    # ----------------------------------------------------------------- #
    
    E_t[p] = 0.5*(I_1*omega_vector[0]**2 + I_2*omega_vector[1]**2 + I_3*omega_vector[2]**2)  # ypologismos ths energeias
    
    # ---------- gia to sfalma ---------- #
    d = E_t[p]-E
    if d == 0:
        sfalma[p] = float("NaN")
       
    else:
        sfalma[p] = m.log10(abs((E_t[p]-E)/E))
    # ----------------------------------- #
    
    dt = dt + h    # allagh xronikou vhmatos gia ton aksona-x twn grafikwn 
    file.write('%f %f %f %f %f\n'%(dt, omega_vector[0], omega_vector[1], omega_vector[2], sfalma[p])) 
   
file.close()

     
data = np.loadtxt('data_project_2_b.txt')
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




