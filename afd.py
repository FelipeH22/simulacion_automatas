import matplotlib.pyplot as plt
from graphviz import Digraph,Source
import numpy as np
import re

class AFD():
    estadosLimbo=list()
    estadosInaccesibles=list()
    Sigma=list()
    Q=list()
    q0=str()
    F=list()
    delta=list()
    def __init__(self,alfabeto,estados,estadoInicial,estadosAceptacion,delta):
        self.Sigma=alfabeto
        self.Q=estados
        self.q0=estadoInicial
        self.F=estadosAceptacion
        self.delta=delta

    def __init__(self,nombreArchivo):
        with open(nombreArchivo) as f:
            contenido=f.readlines()
        contenido=" ".join(list(map(lambda x: x.strip(),contenido))).split("#")
        while '' in contenido: contenido.remove('')
        contenido=list(map(lambda x: x.strip(),contenido))
        for elemento in contenido:
            if '!' in elemento:
                if elemento[1:].strip()!='dfa': raise Exception('El archivo leído no corresponde a un AFD')
            if 'alphabet' in elemento: self.Sigma=sorted(elemento.split(' ')[1:])
            if 'states' in elemento: self.Q=sorted(elemento.split(' ')[1:])
            if 'initial' in elemento: self.q0=elemento.split(' ')[1]
            if 'accepting' in elemento: self.F=sorted(elemento.split(' ')[1:])
            if 'transitions' in elemento: self.delta=sorted(elemento.split(' ')[1:])
        if '-' in self.Sigma[0]: self.Sigma=self.rangoLenguaje(self.Sigma[0])
        self.delta,self.Q=self.creacionDelta(self.delta,self.Sigma,self.Q)
        self.estadosLimbo=self.hallarEstadosLimbo(self.delta,self.Q)
        self.estadosInaccesibles=self.hallarEstadosInaccesibles(self.delta,self.Q,self.q0)

    def rangoLenguaje(self,lenguaje):
        return [chr(caracter) for caracter in range(ord(lenguaje[0]),ord(lenguaje[-1])+1)]

    def manejarTransiciones(self,transiciones):
        return [re.split(':|>', x) for x in transiciones]

    def creacionDelta(self, transiciones, alfabeto, estados):
        matriz=np.asarray([['empt' for _ in range(len(alfabeto))] for _ in range(len(estados))])
        transiciones=self.manejarTransiciones(transiciones)
        for x in transiciones:
            if x[2] not in estados: raise Exception(f"La transición {x[0]}({x[1]}) lleva al estado {x[2]} que no existe entre los estados creados")
            matriz[estados.index(x[0]),alfabeto.index(x[1])]=x[2]
        matriz,estados=self.verificarCorregirCompletitudAFD(matriz,estados)
        return matriz,estados

    def verificarCorregirCompletitudAFD(self,matriz,estados):
        if 'empt' in matriz:
            estados.append('l')
            matriz[matriz=='empt']='l'
            matriz=np.vstack((matriz,['l','l']))
        return matriz,estados

    def hallarEstadosLimbo(self,matriz,estados):
        estados_limbo=list()
        for i in range(len(matriz)):
            if np.array_equal(matriz[i],np.array(['l','l'])): estados_limbo.append(estados[i])
        return estados_limbo

    def hallarEstadosInaccesibles(self,matriz,estados,estado_inicial):
        estados_inaccesibles=list()
        for x in estados:
            if x not in matriz and x!=estado_inicial: estados_inaccesibles.append(x)
        return estados_inaccesibles

    def retornarTransiciones(self,matriz,estados,alfabeto):
        transiciones=list()
        for i in range(len(estados)):
            for j in range(len(alfabeto)):
                transiciones.append(f"{estados[i]}:{alfabeto[j]}>{matriz[i,j]}")
        return transiciones

    def imprimirAFDSimplificado(self):
        inaccesibles=self.hallarEstadosInaccesibles(self.delta,self.Q)
        estados=[x for x in self.Q if x!='l' and x not in inaccesibles]
        transiciones=self.retornarTransiciones(self.delta,self.Q,self.Sigma)
        transiciones=[x for x in transiciones if '>l' not in x]
        print("#!dfa\n"+"#alphabet\n"+"\n".join(self.Sigma)+"\n#states\n"+"\n".join(estados)+"\n#initial\n"+f"{self.q0}\n"+ \
        "#accepting\n"+"\n".join(self.F)+"\n#transitions\n"+"\n".join(transiciones))

    def exportar(self,nombreArchivo):
        contenido=self.toString()
        with open(f"{nombreArchivo}.dfa",'w') as f:
            f.write(contenido)

    def toString(self):
        return "#!dfa\n"+"#alphabet\n"+"\n".join(self.Sigma)+"\n#states\n"+"\n".join(self.Q)+"\n#initial\n"+f"{self.q0}\n"+ \
        "#accepting\n"+"\n".join(self.F)+"\n#transitions\n"+"\n".join(self.retornarTransiciones(self.delta,self.Q,self.Sigma))

    def mostrarGrafo(self,transiciones,estados,estado_inicial,estados_finales,alfabeto):
        transiciones=self.manejarTransiciones(self.retornarTransiciones(transiciones,estados,alfabeto))
        g=Digraph(strict=False)
        g.attr(rankdir="LR")
        g.node('i',label='',shape="point")
        g.edge('i',estado_inicial,arrowsize='0.5')
        for x in estados:
            if x in estados_finales: g.node(x,shape="doublecircle")
            else: g.node(x,shape="circle")
        for x in transiciones:
            g.edge(x[0],x[2],label=x[1],arrowsize='0.5')
        s=Source(g.source,filename="grafo", format="svg")
        s.view()

    def __str__(self):
        self.mostrarGrafo(self.delta,self.Q,self.q0,self.F,self.Sigma)
        return self.toString()



a=AFD('afd.dfa')
print(a)