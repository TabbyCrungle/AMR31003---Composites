import numpy as np

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
s = calc_comp_matrix(E1, E2, G12)
print(E1)
print(s)