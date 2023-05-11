import random

class alfabeto():
    def __init__(self,simbolos):
        self.simbolos=simbolos
    def generarCadenaAleatoria(self,n):
        return " ".join([str(random.choice(self.simbolos)) for _ in range(n)])
