class material:
    def __init__(self, E, V, v) -> None:
        self.E = E
        self.V = V
        self.v = v
        pass

def calc_e1(Ef, Vf, Em, Vm):
    tmp = Ef * Vf + Em * Vm
    return tmp

fibre = material(200000, 0.6, 0.3)
matrix = material(190000, 0.7, 0.1)

E1 = calc_e1(fibre.E,fibre.V,matrix.E,matrix.V) 
print(E1)