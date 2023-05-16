import matplotlib.pyplot as plt
from graphviz import Digraph,Source
from afd import AFD
import numpy as np
import copy
import re

class AFN(AFD):

    def limpiaTransiciones(self,transiciones):
        new_transiciones=list()
        for i in range(len(transiciones)):
            if ';' in transiciones[i][2]:
                division=transiciones[i][2].split(';')
                for x in division:
                    new_transiciones.append([transiciones[i][0],transiciones[i][1],x])
            else: new_transiciones.append(transiciones[i])
        return new_transiciones

    def hallarEstadosInaccesibles(self):
        transiciones=self.manejarTransiciones(self.retornarTransiciones())
        transiciones=self.limpiaTransiciones(transiciones)
        landing_states=[x[2] for x in transiciones]
        return [x for x in self.Q if x not in landing_states]


    def retornarTransiciones(self):
        transiciones=list()
        for i in range(len(self.Q)):
            for j in range(len(self.Sigma)):
                transiciones.append(f"{self.Q[i]}:{self.Sigma[j]}>{self.delta[i,j]}")
        return transiciones

    def imprimirAFDSimplificado(self):
        self.mostrarGrafoSimplificado()
        inaccesibles=self.hallarEstadosInaccesibles()
        estados=[x for x in self.Q if x!='l' and x not in inaccesibles]
        transiciones=self.retornarTransiciones()
        transiciones=[x for x in transiciones if '>l' not in x]
        print("#!nfa\n"+"#alphabet\n"+"\n".join(self.Sigma)+"\n#states\n"+"\n".join(estados)+"\n#initial\n"+f"{self.q0}\n"+ \
        "#accepting\n"+"\n".join(self.F)+"\n#transitions\n"+"\n".join(transiciones))

    def exportar(self,nombreArchivo):
        contenido=self.toString()
        with open(f"{nombreArchivo}.nfa",'w') as f:
            f.write(contenido)

    def toString(self):
        self.mostrarGrafo()
        return "#!nfa\n"+"#alphabet\n"+"\n".join(self.Sigma)+"\n#states\n"+"\n".join(self.Q)+"\n#initial\n"+f"{self.q0}\n"+ \
        "#accepting\n"+"\n".join(self.F)+"\n#transitions\n"+"\n".join(self.retornarTransiciones())

    def mostrarGrafoSimplificado(self,archivo="grafo_s"):
        transiciones=self.manejarTransiciones(self.retornarTransiciones())
        g_s=Digraph(strict=False)
        g_s.attr(rankdir="LR")
        g_s.node('i',label='',shape="point")
        g_s.edge('i',self.q0,arrowsize='0.5')
        for x in self.Q:
            if x in self.F:
                g_s.node(x,shape="doublecircle")
            elif x!='l' and x not in self.estadosInaccesibles:
                g_s.node(x,shape="circle")
        new_transiciones=self.limpiaTransiciones(transiciones)
        for x in new_transiciones:
            if 'l' not in x and x[0] not in self.estadosInaccesibles and x[2] not in self.estadosInaccesibles: g_s.edge(x[0],x[2],label=x[1],arrowsize='0.5')
        s=Source(g_s.source,filename=archivo,format="svg")
        s.view()

    def mostrarGrafo(self,archivo="grafo"):
        transiciones=self.manejarTransiciones(self.retornarTransiciones())
        g=Digraph(strict=False)
        g.attr(rankdir="LR")
        g.node('i',label='',shape="point")
        g.edge('i',self.q0,arrowsize='0.5')
        for x in self.Q:
            if x in self.F: g.node(x,shape="doublecircle")
            else: g.node(x,shape="circle")
        new_transiciones=self.limpiaTransiciones(transiciones)
        for x in new_transiciones:
            g.edge(x[0],x[2],label=x[1],arrowsize='0.5')
        s=Source(g.source,filename=archivo,format="svg")
        s.view()

    def procesarCadena(self,cadena):
        estado_actual = self.q0
        for elemento in cadena:
            estado_actual = self.delta[self.Q.index(estado_actual), self.Sigma.index(elemento)]
        if estado_actual in self.F: return True
        return False

    def procesarCadenaConDetalles(self,cadena):
        estado_actual=self.q0
        cadena_=cadena
        for elemento in cadena:
            if cadena_!='': print(f"[{estado_actual},{cadena_}]", end=" -> ")
            estado_actual=self.delta[self.Q.index(estado_actual),self.Sigma.index(elemento)]
            cadena_=cadena_[1:]
        if estado_actual in self.F:
            print("Aceptación")
            return True
        print("No Aceptación")
        return False

    def procesarListaCadenas(self,listaCadenas,nombreArchivo,imprimirPantalla):
        proceso=list()
        for cadena in listaCadenas:
            resultado=str()
            estado_actual=self.q0
            cadena_=cadena
            for elemento in cadena:
                if cadena_!='':resultado+=f"[{estado_actual},{cadena_}] -> "
                estado_actual=self.delta[self.Q.index(estado_actual),self.Sigma.index(elemento)]
                cadena_=cadena_[1:]
            if estado_actual in self.F:resultado+='si'
            else: resultado+='no'
            proceso.append(resultado)
            if imprimirPantalla: print(resultado)
        try:
            with open(nombreArchivo,'w') as file:
                file.write("\n".join(proceso))
        except:
            with open("salidaProcesamiento.txt",'w') as file:
                file.write("\n".join(proceso))
                print(f"El nombre de archivo '{nombreArchivo}' es invalido, guardando en salidaProcesamiento.txt")

    def hallarComplemento(self,afdInput):
        b=copy.deepcopy(afdInput)
        for estado in b.Q:
            if estado not in b.F: b.F.append(estado)
            else: b.F.remove(estado)
        return b

    def __str__(self): return self.toString()

a=AFN('afn.nfa')
print(a)
print(a.hallarEstadosInaccesibles())
a.mostrarGrafoSimplificado()