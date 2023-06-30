import matplotlib.pyplot as plt
from graphviz import Digraph, Source
from collections import deque
from afd import AFD
import numpy as np
import copy
import re


class AFPD(AFD):
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
        return [re.split(':|>', x) for x in transiciones]

    def creacionDelta(self):
        matriz = np.asarray([['empt' for _ in range(len(self.Sigma))] for _ in range(len(self.Q))],
                            dtype=np.dtype('U100'))
        transiciones = self.manejarTransiciones(self.delta)
        for x in transiciones:
            matriz[self.Q.index(x[0]), self.Sigma.index(x[1])] = f"{x[3]};{x[2]}:{x[4]}"
        self.delta = matriz
        self.verificarCorregirCompletitudAFPD()

    def verificarCorregirCompletitudAFPD(self):
        if 'empt' in self.delta:
            self.Q.append('l')
            self.delta[self.delta == 'empt'] = 'l;$:$'
            self.delta = np.vstack((self.delta, ['l;$:$' for _ in range(len(self.Sigma))]))
        return self.delta, self.Q

    def hallarEstadosLimbo(self):
        estados_limbo = list()
        for i in range(len(self.delta)):
            if np.array_equal(self.delta[i], np.array(['l;$:$' for _ in range(len(self.Sigma))])): estados_limbo.append(
                self.Q[i])
        return estados_limbo

    def hallarEstadosInaccesibles(self):
        estados_accesibles = {self.q0}
        for i in range(len(self.delta)):
            for j in range(len(self.delta[i])):
                estado=self.delta[i,j].split(';')[0]
                if 'l' not in estado and i!=self.Q.index(estado): estados_accesibles.add(estado)
        return list(set(self.Q).difference(estados_accesibles))

    def retornarTransiciones(self):
        transiciones=list()
        for i in range(len(self.Q)):
            for j in range(len(self.Sigma)):
                tr=self.delta[i,j].split(';')
                if 'l' not in tr: transiciones.append([self.Q[i],self.Sigma[j],tr[0],tr[1]])
                else: transiciones.append([self.Q[i],self.Sigma[j],'l',"$:$"])
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
            if 'l' not in x and x[0] not in self.estadosInaccesibles and x[2] not in self.estadosInaccesibles and x[2] not in self.estadosLimbo:
                pila = x[3].split(':')
                pila = [elemento if elemento != '$' else '\u03BB' for elemento in pila]
                g_s.edge(x[0], x[2], label=f"{x[1]},{'|'.join(pila)}", arrowsize='0.5')
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
            g.edge(x[0], x[2], label=f"{x[1]},{'|'.join(x[3].split(':'))}", arrowsize='0.5')
        s = Source(g.source, filename=archivo, format="svg")
        s.view()

    def modificarPila(self, operacion, parametro):
        if operacion == "insercion":
            self.Stack.append(parametro)
        elif operacion == "eliminacion":
            if len(self.Stack) < 1:
                raise ValueError("La pila está vacía")
            elif self.Stack[-1] != parametro:
                raise ValueError("El tope de la pila no corresponde con el parámetro")
            else:
                self.Stack.pop()

    def procesarCadena(self, cadena):
        current_state = self.q0
        while cadena:
            char = cadena[0]
            transicion = self.delta[self.Q.index(current_state), self.Sigma.index(char)]
            current_state, out, inp = re.split(";|:", transicion)
            if out != '$':
                try:
                    self.modificarPila("eliminacion", out)
                except ValueError:
                    break
            if inp != '$': self.modificarPila("insercion", inp)
            cadena = cadena[1:]
        if current_state in self.F and len(self.Stack) == 0 and len(cadena) == 0:
            self.Stack=list()
            return True
        self.Stack=list()
        return False

    def procesarCadenaConDetalles(self, cadena):
        current_state = self.q0
        while cadena:
            print(f"({current_state},{cadena},{''.join(self.Stack) if len(self.Stack) > 0 else '$'})", end="->")
            char = cadena[0]
            transicion = self.delta[self.Q.index(current_state), self.Sigma.index(char)]
            current_state, out, inp = re.split(";|:", transicion)
            if out != '$':
                try:
                    self.modificarPila("eliminacion", out)
                except ValueError:
                    break
            if inp != '$': self.modificarPila("insercion", inp)
            cadena = cadena[1:]
        print(
            f"({current_state},{cadena if len(cadena) > 0 else '$'},{''.join(self.Stack) if len(self.Stack) > 0 else '$'})",
            end="")
        if current_state in self.F and len(self.Stack) == 0 and len(cadena) == 0:
            print(">>accepted")
            self.Stack=list()
            return True
        self.Stack=list()
        print(">>rejected")
        return False

    def procesarListaCadenas(self, listaCadenas, nombreArchivo, imprimirPantalla):
        with open(nombreArchivo,"w") as f:
            for cadena in listaCadenas:
                salida = str()
                current_state = self.q0
                while cadena:
                    salida += f"({current_state},{cadena},{''.join(self.Stack) if len(self.Stack) > 0 else '$'})->"
                    char = cadena[0]
                    transicion = self.delta[self.Q.index(current_state), self.Sigma.index(char)]
                    current_state, out, inp = re.split(";|:", transicion)
                    if out != '$':
                        try:
                            self.modificarPila("eliminacion", out)
                        except ValueError:
                            break
                    if inp != '$': self.modificarPila("insercion", inp)
                    cadena = cadena[1:]
                salida += f"({current_state},{cadena if len(cadena) > 0 else '$'},{''.join(self.Stack) if len(self.Stack) > 0 else '$'})"
                if current_state in self.F and len(self.Stack) == 0 and len(cadena) == 0:
                    salida += ">>accepted"
                else: salida += ">>rejected"
                f.write(salida+"\n")
                if imprimirPantalla: print(salida)
                self.Stack=list()

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
        return "#!dpda\n" + "#states\n" + "\n".join(self.Q) + "\n#initial\n" + f"{self.q0}\n" + \
            "#accepting\n" + "\n".join(self.F) + "\n#tapeAlphabet\n" + "\n".join(self.Sigma) + "#StackAlphabet\n" + "\n".join(self.Gamma) + "\n#transitions\n" + "\n".join(transiciones)

    def __str__(self):
        self.mostrarGrafoSimplificado("_afdp")
        return self.toString()