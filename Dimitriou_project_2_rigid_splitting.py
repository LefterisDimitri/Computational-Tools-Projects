import scipy.special as sp
import numpy as np 
import math as m
import matplotlib.pyplot as plt 

# ---------- STHN ASKHSH AUTH THA LUSOUME TO PROBLHMA TOU STEREOU SWMATOS ME XRHSH THS METHODOU SPLITTING ---------- #


# orismos metablhtwn kai fusikwn statherwn ws global metablhtes

I_1 = 0.8 # roph adraneias ---> I1
I_2 = 0.9 # roph adraneias ---> I2
I_3 = 1.0 # roph adraneias ---> I3

dt = 0.02 # xroniko vhma pou menei stathero se kathe epanalhpsh logo ths methodou pou xrhsimopoihsame
t = 0.02 # xroniko vhma pou mas vohtha stis grafikes parastaseis

#omega_vector = np.array([1.0,0.0,2.0])
M = np.array([0.0,0.0,0.0])    # arxikopoihsh tou voithitikou pinaka M = (M1,M2,M3) 

Mxn = np.array([0.0,0.0,0.0]) # arxikopoihsh twn voithitikwn pinakwn gia tis sunarthseis - HM1: Mxn = (M1,M2,M3) sumfwna me thn theoria
Myn = np.array([0.0,0.0,0.0]) # arxikopoihsh twn voithitikwn pinakwn gia tis sunarthseis - HM2: Myn = (M1,M2,M3) sumfwna me thn theoria
Mzn = np.array([0.0,0.0,0.0]) # arxikopoihsh twn voithitikwn pinakwn gia tis sunarthseis - HM3: Mzn = (M1,M2,M3) sumfwna me thn theoria

Rx = np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]) # arxikopoihsh pinaka strofhs ws pros aksona-x
Ry = np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]) # arxikopoihsh pinaka strofhs ws pros aksona-y
Rz = np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]) # arxikopoihsh pinaka strofhs ws pros aksona-z

result_x = np.array([0.0,0.0,0.0]) # gia ton pollaplasiasmo twn pinakwn
result_y = np.array([0.0,0.0,0.0]) # gia ton pollaplasiasmo twn pinakwn
result_z = np.array([0.0,0.0,0.0]) # gia ton pollaplasiasmo twn pinakwn

omega_n_1 = np.zeros(5000)
omega_n_2 = np.zeros(5000)
omega_n_3 = np.zeros(5000)

E_t = np.zeros(5000)
sfalma = np.zeros(5000)

file = open('data_project_2_c.txt','w')
omega_1_t_0 = 1.0
omega_2_t_0 = 0.0
omega_3_t_0 = 2.0
file.write('%f %f %f %f %f\n'%(0.0, omega_1_t_0, omega_2_t_0, omega_3_t_0, float("NaN")))

E = 0.5*(I_1*omega_1_t_0**2 + I_2*omega_2_t_0**2 + I_3*omega_3_t_0**2)   # energeia gia t = 0 ---> theorhtika statherh

M[0] = 0.8 # M[0] = I1*omega_1 gia t = 0
M[1] = 0.0 # M[1] = I2*omega_2 gia t = 0
M[2] = 2.0 # M[2] = I3*omega_3 gia t = 0

def R_x(M1,Mxn):
    th = (M1[0]/I_1)*dt/2.0
    Rx = ([[1.0,0.0,0.0],[0.0,np.cos(th),np.sin(th)],[0.0,-np.sin(th),np.cos(th)]])
    result_x = np.dot(Rx,M1)
    for i in range(0,3):
        Mxn[i] = result_x[i]
    return Mxn
    #for i in range(0, 3):
        #for j in range(0, 3):
            #Mxn[i] += Rx[i][j]*M[j]
    
def R_y(M1,Myn):
    th = (M1[1]/I_2)*dt/2.0
    Ry = ([[np.cos(th),0.0,-np.sin(th)],[0.0,1.0,0.0],[np.sin(th),0.0,np.cos(th)]])
    result_y = np.dot(Ry,M1)
    for i in range(0,3):
        Myn[i] = result_y[i]
    return Myn
    #for i in range(0, 3):
        #for j in range(0, 3):
            #Myn[i] += Ry[i][j]*M[j]
         
def R_z(M1,Mzn):  
    th = (M1[2]/I_3)*dt
    Rz = ([[np.cos(th),np.sin(th),0.0],[-np.sin(th),np.cos(th),0.0],[0.0,0.0,1.0]])
    result_z = np.dot(Rz,M1)
    for i in range(0,3):
        Mzn[i] = result_z[i]
    return Mzn
    #for i in range(0, 3):
        #for j in range(0, 3):
            #Mzn[i] += Rz[i][j]*M[j]
    
for p in range(1,5000):
    # --- Methodos splitting --- #
    R_x(M,Mxn)     
    R_y(Mxn,Myn)
    R_z(Myn,Mzn)
    R_y(Mzn,Myn)
    R_x(Myn,Mxn)
    M[0] = Mxn[0]                   # ananewsh arxikwn sunthikwn gia thn epomenh epanalhspsh ---> M[0] = M1
    M[1] = Mxn[1]                   # ananewsh arxikwn sunthikwn gia thn epomenh epanalhspsh ---> M[1] = M2
    M[2] = Mxn[2]                   # ananewsh arxikwn sunthikwn gia thn epomenh epanalhspsh ---> M[2] = M3
    # -------------------------- #
    
    # ----- ypologismos (omega_1, omega_2, omega_3) sto telos kathe epanalhpshs, gia kathe dt = 0.02 ----- #
    omega_n_1[p] = Mxn[0]/I_1
    omega_n_2[p] = Mxn[1]/I_2
    omega_n_3[p] = Mxn[2]/I_3
    # ---------------------------------------------------------------------------------------------------- #
    
    E_t[p] = 0.5*(I_1*omega_n_1[p]**2 + I_2*omega_n_2[p]**2 + I_3*omega_n_3[p]**2)   # ypologismos ths energeias
    
    # ---------- gia to sfalma ---------- #
    d = E_t[p]-E
    if d == 0.0:
        sfalma[p] = float("NaN")
       
    else:
        sfalma[p] = m.log10(abs((E_t[p]-E)/E))
    # ----------------------------------- #

    file.write('%f %f %f %f %f\n'%(t, omega_n_1[p], omega_n_2[p], omega_n_3[p], sfalma[p]))
    #file.write('%f %f %f %f %f\n'%(t, M[0], M[1], M[2], sfalma[p]))
    
        
    t = t + 0.02
    Mxn = np.array([0.0,0.0,0.0])
    Myn = np.array([0.0,0.0,0.0])
    Mzn = np.array([0.0,0.0,0.0])
file.close()
    
data = np.loadtxt('data_project_2_c.txt')
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
