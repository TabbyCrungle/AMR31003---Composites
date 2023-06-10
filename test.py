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
    

def calc_e1(Ef, Vf, Em, Vm):
    tmp = Ef * Vf + Em * Vm
    return tmp

fibre = material(200000, 0.6, 0.3)
matrix = material(190000, 0.7, 0.1)

E1 = calc_e1(fibre.E,fibre.V,matrix.E,matrix.V) 
print(E1)