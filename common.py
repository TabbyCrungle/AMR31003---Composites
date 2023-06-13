import numpy as np
from numpy import matlib as mb
from numpy import linalg as lg
import matplotlib.pyplot as plt
import json

  
class material:
    def __init__(self, E, V, v) -> None:
        self.E = E
        self.V = V
        self.v = v
        pass   
    
def calc_g(Ex, vx):
    tmp = (Ex / ( 2 * (1 + vx)))
    return tmp

def calc_e1(Ef, Vf, Em, Vm):
    tmp = (Ef * Vf) + (Em * Vm)
    return tmp

def calc_e2(Ef, Em, Vf, Vm):
    tmp = (Ef * Em) / ((Vf * Em) + (Vm * Vf))
    return tmp 

def calc_g12(Gf, Gm, Vf, Vm):
    tmp = (Gf * Gm) / ((Vf* Gm) + (Vm * Gf))
    return tmp 

def calc_ht(Xf, Xm, ht, Vf):
    eta = ((Xf/ Xm)- 1)/((Xf / Xm) + ht) 
    tmp = (Xm * (1 + (ht * eta * Vf)))/(1 - (eta * Vf))
    return tmp 

def calc_comp_matrix(E1, E2, G12, vf, vm, Vf, Vm):
    s11 = 1/E1
    s22 = 1/E2 
    nu12 = (vf * Vf) + (vm * Vm)
    nu12 = 0.32
    nu21 = (nu12 * E2)/ E1
    s12 = -nu21 / E2
    s66 = 1/G12
    s = np.array([[s11, s12,0],[s12,s22,0],[0,0,s66]])
    return s

jf = "inputs.json"

with open(jf) as data_file:
    inputs = json.load(data_file)
    
try:
    n_layer = int(inputs['parameters']['tlayers'])
except (TypeError, KeyError):
    print("could not read the total number of layers")

fselect = 0
mselect = 1

try:
    inputef = int(inputs['mech props'][fselect]['E'])
    inputvf= inputs['mech props'][fselect]['volfract']
    inputpf = inputs['mech props'][fselect]['poisson']
except (TypeError, KeyError):
    print("could not read fibre properties")    
    
try:
    inputem = int(inputs['mech props'][mselect]['E'])
    inputvm= inputs['mech props'][mselect]['volfract']
    inputpm = inputs['mech props'][mselect]['poisson']
except (TypeError, KeyError):
    print("could not read matrix properties")   
      
fibre = material(inputef, inputvf, inputpf)
matrix = material(inputem, inputvf, inputpf)

try:
    ht = inputs['constants']['halpin-tsai']
except (TypeError, KeyError):
    print("could not read constants") 
    
E1 = calc_e1(fibre.E,fibre.V,matrix.E,matrix.V) 
E2 = calc_e2(fibre.E,matrix.E,fibre.V,matrix.V)
Gf = calc_g(fibre.E,fibre.v)
Gm = calc_g(matrix.E, matrix.v)
G12 = calc_g12(Gf, Gm, fibre.V, matrix.V)
G12 = 7
E2HT = calc_ht(fibre.E, matrix.E, ht, fibre.V)
G12HT = calc_ht(Gf, Gm, ht, fibre.V)

s = np.mat(calc_comp_matrix(E1, E2HT, G12HT, fibre.v, matrix.v, fibre.V, matrix.V))
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
for j in range(0,n_layer):
    try:
        angle = int(inputs['layerlist'][j]['angle'])
    except (TypeError, KeyError):
        print("could not read the angle")
    t_s = s_transformed(s, angle)
    t_q = q_transformed(q, angle)
    q_list.append(t_q)
    
empty = np.zeros((3,3))
AO = np.mat(empty)
BO = np.mat(empty)
DO = np.mat(empty)

try:
    t = float(inputs['parameters']['tthickness']) 
except (TypeError, KeyError):
    print("could not read total thickness")
    
try:
    t_k = int(inputs['layerlist'][0]['thickness'])
except (TypeError, KeyError):
    print("could not read the layer thickness")
z_list = []
tmp4 = 0
tmp5 = -(n_layer)/2*t_k

for i in range(0,(n_layer +1)): 
    try:
        t_k = int(inputs['layerlist'][i-1]['thickness'])
    except (TypeError, KeyError):
        print("could not read the layer thickness")
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


# a = lg.inv(A)

# E1C = 1/a[0,0]
# E2C = 1/a[1,1]
# G12C = 1/a[2,2]
# nu12c = -a[0,1]/a[0,0]
# nu21c = -a[0,1]/a[1,1]
          

print(E1)
print(E2HT)


print("Compliance matrix [S]:")
print(s)
print("Stiffness matrix [Q]:")
print(q)
# print("[A]:")
# print(A)
# print("[B]:")
# print(B)
# print("[D]:")
# print(D)
# print("[a] = [A]^-1:")
# print(a)
# print("Composite longitudinal modulus of elasticity:", E1C)
# print("Composite transverse modulus of elasticity:", E2C)
# print("Composite shear modulus:", G12C)
