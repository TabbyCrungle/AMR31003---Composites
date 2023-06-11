import numpy as np
from numpy import matlib as mb
from numpy import linalg as lg
import matplotlib.pyplot as plt

class material:
    def __init__(self, E, V, v) -> None:
        self.E = E
        self.V = V
        self.v = v
        pass
    
class lamina:
    def __init__(self, E1, E2, G12, E) -> None:
        self.E1 = E1
        self.E2 = E2
        self.G12 = G12
        self.E = E
        pass
    
class layercount:
    def __init__(self, n1, n2, n3) -> None:
        self.n1 = n1
        self.n2 = n2
        self.n3 = n3
        pass
    
class layup:
    def __init__(self, theta1, theta2, theta3) -> None:
        self.theta1 = theta1
        self.theta2 = theta2
        self.theta3 = theta3
        pass
    
class composite:
    def __init__(self, Ec1, Ec2, Gc12) -> None:
        self.Ec1 = Ec1
        self.Ec2 = Ec2
        self.Gc12 = Gc12
        pass    
    
def calc_g(Ex, vx):
    tmp = Ex / 2 * (1 + vx)
    return tmp

def calc_e1(Ef, Vf, Em, Vm):
    tmp = Ef * Vf + Em * Vm
    return tmp

def calc_e2(Ef, Em, Vf, Vm):
    tmp = (Ef * Em) / (Vf * Em + Vm * Vf)
    return tmp 

def calc_g12(Gf, Gm, Vf, Vm):
    tmp = (Gf * Gm) / (Vf* Gm + Vm * Gf)
    return tmp 

def calc_ht(Xf, Xm, ht, Vf, Vm):
    eta = (Xf/ Xm- 1)/(Xf / Xm + ht) 
    tmp = Xm * (1 + (ht * eta * Vf))/1 - (eta * Vf)
    return tmp 

def calc_comp_matrix(E1, E2, G12):
    s11 = 1/E1
    s22 = 1/E2 
    nu12 = (2 * G12)/E1 - 1
    nu21 = (nu12 * E2)/ E1
    s12 = -1 * nu21 / E2
    s66 = 1/G12
    s = np.array([[s11, s12,0],[s12,s22,0],[0,0,s66]])
    return s

fibre = material(200000, 0.6, 0.3)
matrix = material(190000, 0.7, 0.1)
ht = 1

E1 = calc_e1(fibre.E,fibre.V,matrix.E,matrix.V) 
E2 = calc_e2(fibre.E,matrix.E,fibre.V,matrix.V)
Gf = calc_g(fibre.E,fibre.v)
Gm = calc_g(matrix.E, matrix.v)
G12 = calc_g12(Gf, Gm, fibre.V, matrix.V)
E2HT = calc_ht(fibre.E, matrix.E, ht, fibre.V, matrix.V)
G12HT = calc_ht(Gf, Gm, ht, fibre.V, matrix.V)

s = np.mat(calc_comp_matrix(E1, E2, G12))
q = lg.inv(s)
    
def s_transformed(s, angle):
    angle = angle * np.pi/180  
    tmp1 = np.cos(angle)
    tmp2 = np.sin(angle)
    T1 = np.mat([[tmp1**2, tmp2**2, 2*tmp1*tmp2],
                    [tmp2**2, tmp1**2, -2*tmp1*tmp2],
                    [-tmp1*tmp2, tmp1*tmp2, tmp1**2-tmp2**2]])
    T2 = np.mat([[tmp1**2, tmp2**2, tmp1*tmp2],
                    [tmp2**2, tmp1**2, -tmp1*tmp2],
                    [-2*tmp1*tmp2, 2*tmp1*tmp2, tmp1**2-tmp2**2]])
    s_transformed = np.linalg.inv(T2) * s * T1
    return s_transformed

def q_transformed(q, angle):
    angle = angle * np.pi/180  
    tmp1 = np.cos(angle)
    tmp2 = np.sin(angle)
    T1 = np.mat([[tmp1**2, tmp2**2, 2*tmp1*tmp2],
                    [tmp2**2, tmp1**2, -2*tmp1*tmp2],
                    [-tmp1*tmp2, tmp1*tmp2, tmp1**2-tmp2**2]])
    T2 = np.mat([[tmp1**2, tmp2**2, tmp1*tmp2],
                    [tmp2**2, tmp1**2, -tmp1*tmp2],
                    [-2*tmp1*tmp2, 2*tmp1*tmp2, tmp1**2-tmp2**2]])
    q_transformed = np.linalg.inv(T1) * q * T2
    return q_transformed

q_list = []
n_layer = int(8)
for j in range(0,n_layer):
    angle = 15
    t_s = s_transformed(s, angle)
    t_q = q_transformed(q, angle)
    q_list.append(t_q)
    
empty = np.zeros((3,3))
AO = np.mat(empty)
BO = np.mat(empty)
DO = np.mat(empty)

t = float(0.2) 
t_k = t/n_layer           
z_list = []
tmp4 = 0
tmp5 = -(n_layer)/2*t_k

for i in range(0,(n_layer +1)): 
    tmp5 = tmp5 + tmp4
    tmp4 = t_k 
    z_list.append(tmp5)

for i in range(3):             
    for j in range(3):
        for k in range(len(q_list)):
            AO[i,j] = q_list[k][i,j]*(z_list[k+1]-z_list[k]) + AO[i,j]
            BO[i,j] = (1/2) *q_list[k][i,j]*((z_list[k+1])**2-(z_list[k])**2) + BO[i,j]
            DO[i,j] = (1/3) *q_list[k][i,j]*((z_list[k+1])**3 -(z_list[k])**3) + DO[i,j]
            
A = np.around(AO,decimals = 3) 
B = np.around(BO,decimals = 3) 
D = np.around(DO,decimals = 3) 

a = lg.inv(A)

E1C = 1/a[0,0]
E2C = 1/a[1,1]
G12C = 1/a[2,2]
nu12c = -a[0,1]/a[0,0]
nu21c = -a[0,1]/a[1,1]
          
#print(E1)
#print(s)
#print(q_list)
print(a)
print(E1C)
