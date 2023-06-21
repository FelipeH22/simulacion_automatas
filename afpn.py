import matplotlib.pyplot as plt
from graphviz import Digraph, Source
from collections import deque
from afpd import AFPD
import numpy as np
import copy
import re


class AFPN(AFPD):
    def __init__(self, **kwargs):
        self.Stack = list()
        if "nombreArchivo" in kwargs:
            with open(kwargs.get("nombreArchivo")) as f:
                contenido = f.readlines()
            contenido = " ".join(list(map(lambda x: x.strip(), contenido))).split("#")
            while '' in contenido: contenido.remove('')
            contenido = list(map(lambda x: x.strip(), contenido))
            for elemento in contenido:
                if 'tapeAlphabet' in elemento: self.Sigma = sorted(elemento.split(' ')[1:])
                if 'states' in elemento: self.Q = sorted(elemento.split(' ')[1:])
                if 'initial' in elemento: self.q0 = elemento.split(' ')[1]
                if 'accepting' in elemento: self.F = sorted(elemento.split(' ')[1:])
                if 'StackAlphabet' in elemento: self.Gamma = sorted(elemento.split(' ')[1:])
                if 'transitions' in elemento: self.delta = sorted(elemento.split(' ')[1:])
            if '-' in self.Sigma[0]:
                self.Sigma = [*self.Sigma, *self.rangoLenguaje(self.Sigma[0])]
                self.Sigma.pop(0)
            self.Sigma = sorted(self.Sigma)
            if '-' in self.Gamma[0]:
                self.Gamma = [*self.Gamma, *self.rangoLenguaje(self.Gamma[0])]
                self.Gamma.pop(0)
            self.Gamma = sorted(self.Gamma)
            self.creacionDelta()
            self.estadosLimbo = self.hallarEstadosLimbo()
            self.estadosInaccesibles = self.hallarEstadosInaccesibles()
        else:
            self.Sigma = kwargs.get("alfabetoCinta")
            self.Q = kwargs.get("estados")
            self.q0 = kwargs.get("estadoInicial")
            self.F = kwargs.get("estadosAceptacion")
            self.delta = kwargs.get("transiciones")
            self.Gamma = kwargs.get("alfabetoPila")
            if None in (self.Sigma, self.Q, self.q0, self.F, self.Gamma): raise ValueError(
                "Faltan parámetros requeridos o están en formato incorrecto")
            if "deltaEnFormato" not in kwargs: self.creacionDelta()
            self.estadosLimbo = self.hallarEstadosLimbo()
            self.estadosInaccesibles = self.hallarEstadosInaccesibles()

    def manejarTransiciones(self, transiciones):
        return [re.split('>|;', x) for x in transiciones]

    def creacionDelta(self):
        matriz = np.array([['empt' for _ in range(len(self.Sigma))] for _ in range(len(self.Q))],
                            dtype=np.dtype('U100'))
        transiciones = self.manejarTransiciones(self.delta)
        for x in transiciones:
            if x[0].split(':')[1]=='$' and '$' not in self.Sigma:
                self.Sigma.append('$')
                matriz = np.hstack((matriz, np.array(['empt' for _ in range(len(self.Q))]).reshape(len(self.Q),1)))
            if matriz[self.Q.index(x[0].split(':')[0]), self.Sigma.index(x[0].split(':')[1])] =='empt':
                matriz[self.Q.index(x[0].split(':')[0]), self.Sigma.index(x[0].split(':')[1])] = ";".join([f"{estado.split(':')[0]},{x[0].split(':')[-1]}:{estado.split(':')[1]}" for estado in x[1:]])
            else:
                matriz[self.Q.index(x[0].split(':')[0]), self.Sigma.index(x[0].split(':')[1])] += ";"+";".join([f"{estado.split(':')[0]},{x[0].split(':')[-1]}:{estado.split(':')[1]}" for estado in x[1:]])
        self.delta = matriz
        self.verificarCorregirCompletitudAFPN()

    def verificarCorregirCompletitudAFPN(self):
        if 'empt' in self.delta:
            self.Q.append('l')
            self.delta[self.delta == 'empt'] = 'l,$:$'
            self.delta = np.vstack((self.delta, ['l,$:$' for _ in range(len(self.Sigma))]))
        return self.delta, self.Q

    def hallarEstadosLimbo(self):
        estados_limbo = list()
        for i in range(len(self.delta)):
            if np.array_equal(self.delta[i], np.array(['l,$:$' for _ in range(len(self.Sigma))])): estados_limbo.append(
                self.Q[i])
        return estados_limbo

    def hallarEstadosInaccesibles(self):
        estados_accesibles = {self.q0}
        for i in range(len(self.delta)):
            for j in range(len(self.delta[i])):
                estados=self.delta[i,j].split(';')
                estados={x.split(',')[0] for x in estados}
                if 'l' in estados: estados.discard('l')
                if self.Q[i] in estados: estados.discard(self.Q[i])
                estados_accesibles=estados_accesibles.union(estados)
        return list(set(self.Q).difference(estados_accesibles))

    def retornarTransiciones(self):
        transiciones=list()
        for i in range(len(self.Q)):
            for j in range(len(self.Sigma)):
                tr=[self.Q[i],self.Sigma[j]]+[x for x in self.delta[i,j].split(';')]
                transiciones.append(tr)
        return transiciones

    def mostrarGrafoSimplificado(self, archivo="grafo_s"):
        transiciones = self.retornarTransiciones()
        g_s = Digraph(strict=False)
        g_s.attr(rankdir="LR")
        g_s.node('i', label='', shape="point")
        g_s.edge('i', self.q0, arrowsize='0.5')
        for x in self.Q:
            if x in self.F and x not in self.estadosInaccesibles:
                g_s.node(x, shape="doublecircle")
            elif x != 'l' and x not in self.estadosInaccesibles and x not in self.estadosLimbo:
                g_s.node(x, shape="circle")
        for x in transiciones:
            for posibles in x[2:]:
                ax=posibles.split(',')
                if x[0] not in self.estadosInaccesibles and ax[0] not in self.estadosInaccesibles and ax[0] not in self.estadosLimbo:
                    pila=ax[1].split(':')
                    pila = [elemento if elemento!='$' else '\u03BB' for elemento in pila]
                    if x[1]=='$': x[1]="\u03BB"
                    g_s.edge(x[0], ax[0], label=f"{x[1]},{'|'.join(pila)}", arrowsize='0.5')
        s = Source(g_s.source, filename=archivo, format="svg")
        s.view()

    def mostrarGrafo(self, archivo="grafo"):
        transiciones = self.retornarTransiciones()
        g = Digraph(strict=False)
        g.attr(rankdir="LR")
        g.node('i', label='', shape="point")
        g.edge('i', self.q0, arrowsize='0.5')
        for x in self.Q:
            if x in self.F:
                g.node(x, shape="doublecircle")
            else:
                g.node(x, shape="circle")
        for x in transiciones:
            for posibles in x[2:]:
                ax = posibles.split(',')
                pila = ax[1].split(':')
                pila = [elemento if elemento != '$' else '\u03BB' for elemento in pila]
                if x[1] == '$': x[1] = "\u03BB"
                g.edge(x[0], ax[0], label=f"{x[1]},{'|'.join(pila)}", arrowsize='0.5')
        s = Source(g.source, filename=archivo, format="svg")
        s.view()

    def procesarCadena(self, cadena):
        if self.procesamiento(cadena)[0]: return True
        return False

    def procesarCadenaConDetalles(self, cadena):
        resultado = self.procesamiento(cadena)
        if resultado[0]: print("->".join(resultado[0][0]),"->Aceptada")
        for i in range(len(resultado[1])):
            print("->".join(resultado[1][i]),"->Rechazada")

    def computarTodosLosProcesamientos(self,cadena,nombreArchivo):
        resultados=self.procesamiento(cadena)
        with open(f"{nombreArchivo}AceptadasAFPN.txt", "w") as file:
            for route in resultados[0]:
                resultado="->".join(route) + "->Aceptación" +"\n"
                file.write(resultado)
                print(resultado)
        with open(f"{nombreArchivo}RechazadasAFPN.txt", "w") as file:
            for route in resultados[1]:
                resultado = "->".join(route) + "->No Aceptación" + "\n"
                file.write(resultado)
                print(resultado)
        with open(f"{nombreArchivo}AbortadasAFPN.txt", "w") as file:
            for route in resultados[2]:
                resultado = "->".join(route) + "->Abortado" + "\n"
                file.write(resultado)
                print(resultado)
        return len(resultado[0])+len(resultado[1])+len(resultado[2])

    def procesarListaCadenas(self, listaCadenas, nombreArchivo, imprimirPantalla):
        with open(nombreArchivo,"w") as f:
            for x in listaCadenas:
                f.write(x+"\n")
                resultado=self.procesarCadenaConDetalles(x)
                if resultado[0]: f.write(resultado[0]+"\n")
                elif resultado[1]: f.write(resultado[1]+"\n")
                else: f.write(resultado[2]+"\n")
                f.write(str(len(resultado[0])+len(resultado[1])+len(resultado[2]))+"\n")
                f.write(str(len(resultado[0]))+"\n")
                f.write(str(len(resultado[1]))+"\n")
                if resultado[0]: f.write("yes")
                else: f.write("no")

    def procesamiento(self, cadena):
        start_state = self.q0
        accepted = list()
        rejected = list()
        aborted = list()
        def dfs(current_states, remaining_input, current_route, stack):
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
                    transitions = self.delta[self.Q.index(next_state), char_index].split(";")
                    for transition in transitions:
                        transition_elements = transition.split(",")
                        possible_state = transition_elements[0]
                        element_to_pop = transition_elements[1].split(":")[0]
                        element_to_add = transition_elements[1].split(":")[1]
                        if element_to_pop != '$':
                            try:
                                self.modificarPila("eliminacion", element_to_pop)
                            except ValueError as e:
                                aborted.append(current_route + [f"['stack_error',{e}]"])
                                continue
                        if element_to_add != '$':
                            try:
                                self.modificarPila("insercion", element_to_add)
                            except ValueError as e:
                                aborted.append(current_route + [f"['stack_error',{e}]"])
                                continue
                        next_route = f"[{possible_state},{remaining_input[:]}]"
                        dfs(possible_state, remaining_input[1:], current_route + [next_route], stack)
                        if element_to_add != '$':
                            self.modificarPila("eliminacion", element_to_add)
                        if element_to_pop != '$':
                            self.modificarPila("insercion", element_to_pop)

        current_route = list()
        dfs(str(start_state), cadena, current_route, [])
        return accepted, rejected, aborted

    def hallarProductoCartesianoConAFD(self, afd):
        aceptacion=list()
        estados=list()
        aux=list()
        transiciones=np.array([['empty' for _ in range(len(self.Sigma))] for _ in range(len(set(self.Q).difference({'l'}))*len(set(afd.Q).difference({'l'})))], dtype=np.dtype('U100'))
        for i in range(len(self.Q)):
            for j in range(len(afd.Q)):
                if 'l' not in self.Q[i] and 'l' not in afd.Q[j]: estados.append(f"{{{self.Q[i]},{afd.Q[j]}}}")
                if self.Q[i] in self.F and afd.Q[j] in afd.F: aceptacion.append(f"{{{self.Q[i]},{afd.Q[j]}}}")
                for z in range(len(self.Sigma)):
                    transicion1=self.delta[i,z].split(';')
                    transicion2=afd.delta[j,z]
                    if 'l' in transicion1 and 'l' in transicion2: new_transicion='l;$:$'
                    elif 'l' in transicion1:
                        new_transicion=f"{transicion2};$:$"
                        if transicion2 not in aux: aux.append(transicion2)
                    elif 'l' in transicion2:
                        new_transicion=f"{transicion1[0]};{transicion1.split(':')[0]}:{transicion1.split(':')[1]}"
                        if transicion1[0] not in aux: aux.append(transicion1[0])
                    else: new_transicion=f"{{{transicion1[0]},{transicion2}}};{transicion1[1]}"
                    transiciones[len(estados)-1,z]=new_transicion
        estados=estados+aux
        for _ in aux: transiciones = np.vstack((transiciones, ['l;$:$' for _ in range(len(self.Sigma))]))
        result=copy.deepcopy(self)
        result.Q=estados
        result.delta=transiciones
        result.q0=f"{{{self.q0},{afd.q0}}}"
        result.F = aceptacion
        result.estadosLimbo=result.hallarEstadosLimbo()
        result.estadosInaccesibles=result.hallarEstadosInaccesibles()
        return result

    def toString(self):
        transiciones=self.retornarTransiciones()
        transiciones=[f"{transicion[0]}:{transicion[1]}:{transicion[3].split(':')[0]}>{transicion[2]}:{transicion[3].split(':')[1]}" for transicion in transiciones]
        return "#!pda\n" + "#states\n" + "\n".join(self.Q) + "\n#initial\n" + f"{self.q0}\n" + \
            "#accepting\n" + "\n".join(self.F) + "\n#tapeAlphabet\n" + "\n".join(self.Sigma) + "#StackAlphabet\n" + "\n".join(self.Gamma) + "\n#transitions\n" + "\n".join(transiciones)

    def __str__(self):
        self.mostrarGrafoSimplificado("_afpn")
        return self.toString()


apn=AFPN(nombreArchivo="./automatas/afpn.pda")
#apn.mostrarGrafoSimplificado()
print(apn.procesarCadena("aaabbb"))
apn.procesarCadenaConDetalles("aaabbb")
apn.computarTodosLosProcesamientos("aaabbb","zzz")

