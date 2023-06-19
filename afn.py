import matplotlib.pyplot as plt
from graphviz import Digraph,Source
from collections import deque
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
        return [x for x in self.Q if x not in landing_states and x!=self.q0]

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
            if x in self.F and x not in self.estadosInaccesibles:
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
        afd = copy.deepcopy(afn)
        alfabeto = afd.Sigma
        estados = afd.Q
        estadoInicial = afd.q0
        estadosAceptacion = afd.F
        transiciones = afd.delta
        previous_size = transiciones.shape
        current_size = (0, 0)
        while previous_size != current_size:
            previous_size = transiciones.shape
            for i in range(len(estados)):
                for j in range(len(alfabeto)):
                    if ";" in transiciones[i, j] and ';'.join(sorted(transiciones[i, j].split(';'))) not in estados:
                        estados.append(';'.join(sorted(transiciones[i, j].split(';'))))
                        transiciones = np.vstack((transiciones, ["empt" for _ in range(len(alfabeto))]))
                        est = set(sorted(transiciones[i, j].split(";")))
                        for x in est:
                            if x in estadosAceptacion: estadosAceptacion.append(';'.join(sorted(transiciones[i, j].split(';'))))
                            for z in range(len(alfabeto)):
                                if transiciones[-1, z] == "empt" or transiciones[-1, z] == 'l':
                                    transiciones[-1, z] = transiciones[estados.index(x), z]
                                elif 'l' not in transiciones[-1, z] and 'l' not in transiciones[estados.index(x), z]:
                                    transiciones[-1, z] = ';'.join(set(transiciones[-1, z].split(';')).union(set(transiciones[estados.index(x), z].split(';'))))
            current_size = transiciones.shape
        afd.delta = transiciones
        for i in range(len(afd.delta)):
            for j in range(len(afd.delta[0])):
                if ';' in afd.delta[i,j]: afd.delta[i,j]=';'.join(sorted(afd.delta[i,j].split(';')))
        afd_equivalente = AFD(alfabeto=alfabeto, estados=estados, estadoInicial=estadoInicial,
                              estadosAceptacion=estadosAceptacion, transiciones=afd.delta, deltaEnFormato=True)
        afd_equivalente.eliminarEstadosInaccesibles()
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

    def procesarCadenaConDetalles(self, cadena):
        resultado = self.procesamientoConDetalles(cadena, self.q0, [])
        if resultado[0]:
            print("->".join(resultado[1]) + "->Aceptación")
        else:
            print("No aceptación")

    def procesamientoConDetalles(self, cadena, current_state, procesamiento):
        if cadena == '':
            return current_state in self.F, procesamiento + [f"[{current_state},]"]

        transiciones = self.delta[self.Q.index(current_state), self.Sigma.index(cadena[0])]
        next_states = [state for state in transiciones.split(';')]

        for next_state in next_states:
            new_procesamiento = procesamiento + [f"[{current_state},{cadena}]"]
            result, path = self.procesamientoConDetalles(cadena[1:], next_state, new_procesamiento)

            if result:
                return True, path

        return False, procesamiento

    def computarTodosLosProcesamientos(self,cadena,nombreArchivo):
        start_state = self.q0
        accepted = list()
        rejected = list()
        aborted = list()
        def dfs(current_states, remaining_input, current_route):
            if remaining_input == "":
                final_states = set(current_states.split(";"))
                if final_states.intersection(self.F):
                    accepted.append(current_route)
                else:
                    rejected.append(current_route)
                return
            current_char = remaining_input[0]
            char_index = self.Sigma.index(current_char)
            for next_state in current_states.split(";"):
                if next_state == 'l':
                    aborted.append(current_route+[f"['l',{remaining_input[:]}]"])
                else:
                    next_states = self.delta[self.Q.index(next_state), char_index]
                    next_route = f"[{next_state},{remaining_input[:]}]"
                    dfs(next_states, remaining_input[1:], current_route + [next_route])
        current_route = list()
        dfs(str(start_state), cadena, current_route)
        with open(f"{nombreArchivo}Aceptadas.txt", "w") as file:
            for route in accepted:
                resultado="->".join(route) + "->Aceptación" +"\n"
                file.write(resultado)
                print(resultado)
        with open(f"{nombreArchivo}Rechazadas.txt", "w") as file:
            for route in rejected:
                resultado = "->".join(route) + "->No Aceptación" + "\n"
                file.write(resultado)
                print(resultado)
        with open(f"{nombreArchivo}Abortadas.txt", "w") as file:
            for route in aborted:
                resultado = "->".join(route) + "->Abortado" + "\n"
                file.write(resultado)
                print(resultado)
        return len(accepted)+len(rejected)+len(aborted)

    def procesarListaCadenas(self,listaCadenas,nombreArchivo,imprimirPantalla):
        start_state = self.q0
        accepted = list()
        rejected = list()
        aborted = list()
        def dfs(current_states, remaining_input, current_route):
            if remaining_input == "":
                final_states = set(current_states.split(";"))
                if final_states.intersection(self.F):
                    accepted.append(current_route)
                else:
                    rejected.append(current_route)
                return
            current_char = remaining_input[0]
            char_index = self.Sigma.index(current_char)
            for next_state in current_states.split(";"):
                if next_state == 'l':
                    aborted.append(current_route + [f"['l',{remaining_input[:]}]"])
                else:
                    next_states = self.delta[self.Q.index(next_state), char_index]
                    next_route = f"[{next_state},{remaining_input[:]}]"
                    dfs(next_states, remaining_input[1:], current_route + [next_route])

        for cadena in listaCadenas:
            current_route = list()
            dfs(str(start_state), cadena, current_route)
            try:
                with open(f"{nombreArchivo}.txt",'a') as file:
                    file.write(cadena+"\n")
                    for route in accepted:
                        resultado = "->".join(route) + "->Aceptación" + "\n"
                        file.write(resultado)
                        if imprimirPantalla: print(resultado)
                    for route in rejected:
                        resultado = "->".join(route) + "->No Aceptación" + "\n"
                        file.write(resultado)
                        if imprimirPantalla: print(resultado)
                    for route in aborted:
                        resultado = "->".join(route) + "->Abortado" + "\n"
                        file.write(resultado)
                        if imprimirPantalla: print(resultado)
            except:
                print(f"El nombre de archivo '{nombreArchivo}' es invalido, guardando en salidaProcesamiento.txt")
                with open(f"{nombreArchivo}.txt", 'a') as file:
                    file.write(cadena + "\n")
                    for route in accepted:
                        resultado = "->".join(route) + "->Aceptación" + "\n"
                        file.write(resultado)
                        if imprimirPantalla: print(resultado)
                    for route in rejected:
                        resultado = "->".join(route) + "->No Aceptación" + "\n"
                        file.write(resultado)
                        if imprimirPantalla: print(resultado)
                    for route in aborted:
                        resultado = "->".join(route) + "->Abortado" + "\n"
                        file.write(resultado)
                        if imprimirPantalla: print(resultado)

    def procesarCadenaConversion(self,cadena):
        a=copy.deepcopy(self)
        b=self.AFNtoAFD(a)
        return b.procesarCadena(cadena)

    def procesarCadenaConDetallesConversion(self,cadena):
        a = copy.deepcopy(self)
        b = self.AFNtoAFD(a)
        return b.procesarCadenaConDetalles(cadena)

    def procesarListaCadenasConversion(self,listaCadenas,nombreArchivo,imprimitPantalla):
        a = copy.deepcopy(self)
        b = self.AFNtoAFD(a)
        return b.procesarListaCadenas(listaCadenas,nombreArchivo,imprimitPantalla)

    def __str__(self):
        self.mostrarGrafoSimplificado("afn")
        return "#!nfa"+self.toString()[5:]
