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

    def imprimirAFNSimplificado(self):
        self.mostrarGrafoSimplificado()
        inaccesibles=self.hallarEstadosInaccesibles()
        estados=[x for x in self.Q if x!='l' and x not in inaccesibles]
        transiciones=self.retornarTransiciones()
        transiciones=[x for x in transiciones if '>l' not in x]
        print("#!nfa\n"+"#alphabet\n"+"\n".join(self.Sigma)+"\n#states\n"+"\n".join(estados)+"\n#initial\n"+f"{self.q0}\n"+ \
        "#accepting\n"+"\n".join(self.F)+"\n#transitions\n"+"\n".join(transiciones))

    def exportar(self,nombreArchivo):
        contenido="#!nfa"+self.toString()[5:]
        with open(f"{nombreArchivo}.nfa",'w') as f:
            f.write(contenido)

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

    def AFNtoAFD(self, afn):
        afd=copy.deepcopy(afn)
        alfabeto=afd.Sigma
        estados=afd.Q
        estadoInicial=afd.q0
        estadosAceptacion=afd.F
        transiciones=afd.delta
        previous_size=transiciones.shape
        current_size=(0,0)
        while previous_size!=current_size:
            previous_size = transiciones.shape
            for i in range(len(estados)):
                for j in range(len(alfabeto)):
                    if ";" in transiciones[i,j] and transiciones[i,j] not in estados:
                        estados.append(transiciones[i,j])
                        transiciones = np.vstack((transiciones, ["empt" for _ in range(len(alfabeto))]))
                        est=transiciones[i,j].split(";")
                        for x in est:
                            for z in range(len(alfabeto)):
                                if transiciones[-1,z]=="empt" or transiciones[-1,z]=='l': transiciones[-1,z]=transiciones[estados.index(x),z]
                                elif transiciones[-1,z]!='l':
                                    if 'l' in transiciones[estados.index(x), z]: continue
                                    transiciones[-1,z]=transiciones[-1,z]+";"+transiciones[estados.index(x),z]
            current_size=transiciones.shape
        afd.delta=transiciones
        transiciones=afd.retornarTransiciones()
        afd_equivalente=AFD(alfabeto=alfabeto,estados=estados,estadoInicial=estadoInicial,estadosAceptacion=estadosAceptacion,transiciones=transiciones)
        afd_equivalente.eliminarEstadosInaccesibles()
        afd_equivalente.eliminarEstadosLimbo()
        return afd_equivalente

    def procesarCadena(self,cadena):
        return self.procesamiento(cadena,self.q0)

    def procesamiento(self,cadena,current_state):
        if cadena == '': return current_state in self.F
        transiciones = self.delta[self.Q.index(current_state), self.Sigma.index(cadena[0])]
        next_states = [state for state in transiciones.split(';')]
        for next_state in next_states:
            if self.procesamiento(cadena[1:],next_state):
                return True
        return False

    def procesarCadenaConDetalles(self,cadena):
        resultado=self.procesamientoConDetalles(cadena,self.q0,list())
        if resultado[0]: print("->".join(reversed(resultado[1]))+"->Aceptación")
        else: print("No aceptación")

    def procesamientoConDetalles(self,cadena,current_state,procesamiento):
        if cadena == '': return current_state in self.F,procesamiento
        transiciones = self.delta[self.Q.index(current_state), self.Sigma.index(cadena[0])]
        next_states = [state for state in transiciones.split(';')]
        for next_state in next_states:
            if self.procesamientoConDetalles(cadena[1:],next_state,procesamiento)[0]:
                procesamiento.append(f"[{next_state},{cadena}]")
                return True,procesamiento
        return False,procesamiento

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

    def __str__(self): return "#!nfa"+self.toString()[5:]

a=AFN(nombreArchivo="afn.nfa")
a.procesarCadenaConDetalles("aaaa")