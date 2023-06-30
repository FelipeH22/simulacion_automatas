import matplotlib.pyplot as plt
import numpy as np
from graphviz import Digraph, Source
import copy
import re


class AFD:
    estadosLimbo = list()
    estadosInaccesibles = list()
    Sigma = list()
    Q = list()
    q0 = str()
    F = list()
    delta = list()

    def __init__(self, **kwargs):
        if "nombreArchivo" in kwargs:
            with open(kwargs.get("nombreArchivo")) as f:
                contenido = f.readlines()
            contenido = " ".join(list(map(lambda x: x.strip(), contenido))).split("#")
            while '' in contenido: contenido.remove('')
            contenido = list(map(lambda x: x.strip(), contenido))
            for elemento in contenido:
                if 'alphabet' in elemento: self.Sigma = sorted(elemento.split(' ')[1:])
                if 'states' in elemento: self.Q = sorted(elemento.split(' ')[1:])
                if 'initial' in elemento: self.q0 = elemento.split(' ')[1]
                if 'accepting' in elemento: self.F = sorted(elemento.split(' ')[1:])
                if 'transitions' in elemento: self.delta = sorted(elemento.split(' ')[1:])
            if '-' in self.Sigma[0]:
                self.Sigma = [*self.Sigma, *self.rangoLenguaje(self.Sigma[0])]
                self.Sigma.pop(0)
            self.Sigma = sorted(self.Sigma)
            if "esLambda" in kwargs:
                if kwargs.get("esLambda") == True: self.Sigma.append('$')
            self.creacionDelta()
            self.estadosLimbo = self.hallarEstadosLimbo()
            self.estadosInaccesibles = self.hallarEstadosInaccesibles()
        else:
            self.Sigma = kwargs.get("alfabeto")
            self.Q = kwargs.get("estados")
            self.q0 = kwargs.get("estadoInicial")
            self.F = kwargs.get("estadosAceptacion")
            self.delta = kwargs.get("transiciones")
            if None in (self.Sigma,self.Q,self.q0,self.F): raise ValueError("Faltan parámetros requeridos o están en formato incorrecto")
            if "deltaEnFormato" not in kwargs: self.creacionDelta()
            self.estadosLimbo = self.hallarEstadosLimbo()
            self.estadosInaccesibles = self.hallarEstadosInaccesibles()


    def rangoLenguaje(self, lenguaje):
        return [chr(caracter) for caracter in range(ord(lenguaje[0]), ord(lenguaje[-1]) + 1)]

    def manejarTransiciones(self, transiciones):
        return [re.split(':|>', x) for x in transiciones]

    def creacionDelta(self):
        matriz = np.asarray([['empt' for _ in range(len(self.Sigma))] for _ in range(len(self.Q))],
                            dtype=np.dtype('U100'))
        transiciones = self.manejarTransiciones(self.delta)
        for x in transiciones:
            matriz[self.Q.index(x[0]), self.Sigma.index(x[1])] = x[2]
        self.delta = matriz
        self.verificarCorregirCompletitudAFD()

    def verificarCorregirCompletitudAFD(self):
        if 'empt' in self.delta:
            self.Q.append('l')
            self.delta[self.delta == 'empt'] = 'l'
            self.delta = np.vstack((self.delta, ['l' for _ in range(len(self.Sigma))]))
        return self.delta, self.Q

    def hallarEstadosLimbo(self):
        estados_limbo = list()
        for i in range(len(self.delta)):
            if np.array_equal(self.delta[i], np.array(['l' for _ in range(len(self.Sigma))])): estados_limbo.append(
                self.Q[i])
        return estados_limbo

    def hallarEstadosInaccesibles(self):
        estados_inaccesibles = list()
        for x in self.Q:
            if (x not in self.delta or (
                    len(np.where(self.delta == x)[0].tolist()) == 1 and np.where(self.delta == x)[0].tolist() == [
                self.Q.index(x)])) and x != self.q0: estados_inaccesibles.append(x)
        return estados_inaccesibles

    def retornarTransiciones(self):
        transiciones = list()
        for i in range(len(self.Q)):
            for j in range(len(self.Sigma)):
                transiciones.append(f"{self.Q[i]}:{self.Sigma[j]}>{self.delta[i, j]}")
        return transiciones

    def imprimirAFDSimplificado(self):
        self.mostrarGrafoSimplificado()
        inaccesibles = self.hallarEstadosInaccesibles()
        estados = [x for x in self.Q if x != 'l' and x not in inaccesibles]
        transiciones = self.retornarTransiciones()
        transiciones = [x for x in transiciones if '>l' not in x]
        print("#!dfa\n" + "#alphabet\n" + "\n".join(self.Sigma) + "\n#states\n" + "\n".join(
            estados) + "\n#initial\n" + f"{self.q0}\n" + \
              "#accepting\n" + "\n".join(self.F) + "\n#transitions\n" + "\n".join(transiciones))

    def exportar(self, nombreArchivo):
        contenido = self.toString()
        with open(f"{nombreArchivo}.dfa", 'w') as f:
            f.write(contenido)

    def toString(self):

        return "#!dfa\n" + "#alphabet\n" + "\n".join(self.Sigma) + "\n#states\n" + "\n".join(
            self.Q) + "\n#initial\n" + f"{self.q0}\n" + \
            "#accepting\n" + "\n".join(self.F) + "\n#transitions\n" + "\n".join(self.retornarTransiciones())

    def mostrarGrafoSimplificado(self, archivo="grafo_s"):
        transiciones = self.manejarTransiciones(self.retornarTransiciones())
        g_s = Digraph(strict=False)
        g_s.attr(rankdir="LR")
        g_s.node('i', label='', shape="point")
        g_s.edge('i', self.q0, arrowsize='0.5')
        for x in self.Q:
            if x in self.F and x not in self.estadosInaccesibles:
                g_s.node(x, shape="doublecircle")
            elif x != 'l' and x not in self.estadosInaccesibles:
                g_s.node(x, shape="circle")
        for x in transiciones:
            if 'l' not in x and x[0] not in self.estadosInaccesibles and x[2] not in self.estadosInaccesibles: g_s.edge(
                x[0], x[2], label=x[1], arrowsize='0.5')
        s = Source(g_s.source, filename=archivo, format="svg")
        s.view()

    def mostrarGrafo(self, archivo="grafo"):
        transiciones = self.manejarTransiciones(self.retornarTransiciones())
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
            g.edge(x[0], x[2], label=x[1], arrowsize='0.5')
        s = Source(g.source, filename=archivo, format="svg")
        s.view()

    def procesarCadena(self, cadena):
        estado_actual = self.q0
        for elemento in cadena:
            estado_actual = self.delta[self.Q.index(estado_actual), self.Sigma.index(elemento)]
        if estado_actual in self.F: return True
        return False

    def procesarCadenaConDetalles(self, cadena):
        estado_actual = self.q0
        cadena_ = cadena
        for elemento in cadena:
            if cadena_ != '': print(f"[{estado_actual},{cadena_}]", end=" -> ")
            estado_actual = self.delta[self.Q.index(estado_actual), self.Sigma.index(elemento)]
            cadena_ = cadena_[1:]
        if estado_actual in self.F:
            print("Aceptación")
            return True
        print("No Aceptación")
        return False

    def procesarListaCadenas(self, listaCadenas, nombreArchivo, imprimirPantalla):
        proceso = list()
        for cadena in listaCadenas:
            resultado = str()
            estado_actual = self.q0
            cadena_ = cadena
            for elemento in cadena:
                if cadena_ != '': resultado += f"[{estado_actual},{cadena_}] -> "
                estado_actual = self.delta[self.Q.index(estado_actual), self.Sigma.index(elemento)]
                cadena_ = cadena_[1:]
            if estado_actual in self.F:
                resultado += 'si'
            else:
                resultado += 'no'
            proceso.append(resultado)
            if imprimirPantalla: print(resultado)
        try:
            with open(f"{nombreArchivo}.txt", 'w') as file:
                file.write("\n".join(proceso))
        except:
            with open("salidaProcesamiento.txt", 'w') as file:
                file.write("\n".join(proceso))
                print(f"El nombre de archivo '{nombreArchivo}' es invalido, guardando en salidaProcesamiento.txt")

    def hallarComplemento(self, afdInput):
        b = copy.deepcopy(afdInput)
        for estado in b.Q:
            if estado not in b.F:
                b.F.append(estado)
            else:
                b.F.remove(estado)
        return b

    def hallarProductoCartesianoY(self, afd1, afd2):
        estados = list()
        estados_finales = list()
        nuevas_transiciones = list()
        resultado = copy.deepcopy(afd1)
        for estado_a1 in afd1.Q:
            for estado_a2 in afd2.Q:
                for elemento in resultado.Sigma:
                    nuevas_transiciones.append(
                        f"{{{estado_a1},{estado_a2}}}:{elemento}>{{{afd1.delta[afd1.Q.index(estado_a1), afd1.Sigma.index(elemento)]},"
                        f"{afd2.delta[afd2.Q.index(estado_a2), afd2.Sigma.index(elemento)]}}}")
                    print(
                        f"\u03B4(({estado_a1},{estado_a2}),{elemento}) = (\u03B4\N{SUBSCRIPT ONE}({estado_a1},{elemento}),\u03B4\N{SUBSCRIPT TWO}({estado_a2},{elemento}))="
                        f"({afd1.delta[afd1.Q.index(estado_a1), afd1.Sigma.index(elemento)]},{afd2.delta[afd2.Q.index(estado_a2), afd2.Sigma.index(elemento)]}),")
                estados.append(f"{{{estado_a1},{estado_a2}}}")
                if estado_a1 in afd1.F and estado_a2 in afd2.F: estados_finales.append(f"{{{estado_a1},{estado_a2}}}")
        resultado.Q = estados
        resultado.F = estados_finales
        resultado.delta = nuevas_transiciones
        resultado.creacionDelta()
        resultado.q0 = f"{{{afd1.q0},{afd2.q0}}}"
        resultado.estadosInaccesibles = resultado.hallarEstadosInaccesibles()
        while resultado.estadosInaccesibles:
            x = resultado.estadosInaccesibles[0]
            resultado.delta = np.delete(resultado.delta, resultado.Q.index(x), 0)
            resultado.Q.remove(x)
            resultado.estadosInaccesibles = resultado.hallarEstadosInaccesibles()
        return resultado

    def hallarProductoCartesianoO(self, afd1, afd2):
        estados = list()
        estados_finales = list()
        nuevas_transiciones = list()
        resultado = copy.deepcopy(afd1)
        for estado_a1 in afd1.Q:
            for estado_a2 in afd2.Q:
                for elemento in resultado.Sigma:
                    nuevas_transiciones.append(
                        f"{{{estado_a1},{estado_a2}}}:{elemento}>{{{afd1.delta[afd1.Q.index(estado_a1), afd1.Sigma.index(elemento)]},"
                        f"{afd2.delta[afd2.Q.index(estado_a2), afd2.Sigma.index(elemento)]}}}")
                    print(
                        f"\u03B4(({estado_a1},{estado_a2}),{elemento}) = (\u03B4\N{SUBSCRIPT ONE}({estado_a1},{elemento}),\u03B4\N{SUBSCRIPT TWO}({estado_a2},{elemento}))="
                        f"({afd1.delta[afd1.Q.index(estado_a1), afd1.Sigma.index(elemento)]},{afd2.delta[afd2.Q.index(estado_a2), afd2.Sigma.index(elemento)]}),")
                estados.append(f"{{{estado_a1},{estado_a2}}}")
                if estado_a1 in afd1.F or estado_a2 in afd2.F: estados_finales.append(f"{{{estado_a1},{estado_a2}}}")
        resultado.Q = estados
        resultado.F = estados_finales
        resultado.delta = nuevas_transiciones
        resultado.creacionDelta()
        resultado.q0 = f"{{{afd1.q0},{afd2.q0}}}"
        resultado.estadosInaccesibles = resultado.hallarEstadosInaccesibles()
        while resultado.estadosInaccesibles:
            x = resultado.estadosInaccesibles[0]
            resultado.delta = np.delete(resultado.delta, resultado.Q.index(x), 0)
            resultado.Q.remove(x)
            resultado.estadosInaccesibles = resultado.hallarEstadosInaccesibles()
        return resultado

    def hallarProductoCartesianoDiferencia(self, afd1, afd2):
        estados = list()
        estados_finales = list()
        nuevas_transiciones = list()
        resultado = copy.deepcopy(afd1)
        for estado_a1 in afd1.Q:
            for estado_a2 in afd2.Q:
                for elemento in resultado.Sigma:
                    nuevas_transiciones.append(
                        f"{{{estado_a1},{estado_a2}}}:{elemento}>{{{afd1.delta[afd1.Q.index(estado_a1), afd1.Sigma.index(elemento)]},"
                        f"{afd2.delta[afd2.Q.index(estado_a2), afd2.Sigma.index(elemento)]}}}")
                    print(
                        f"\u03B4(({estado_a1},{estado_a2}),{elemento}) = (\u03B4\N{SUBSCRIPT ONE}({estado_a1},{elemento}),\u03B4\N{SUBSCRIPT TWO}({estado_a2},{elemento}))="
                        f"({afd1.delta[afd1.Q.index(estado_a1), afd1.Sigma.index(elemento)]},{afd2.delta[afd2.Q.index(estado_a2), afd2.Sigma.index(elemento)]}),")
                estados.append(f"{{{estado_a1},{estado_a2}}}")
                if estado_a1 in afd1.F and estado_a2 not in afd2.F: estados_finales.append(
                    f"{{{estado_a1},{estado_a2}}}")
        resultado.Q = estados
        resultado.F = estados_finales
        resultado.delta = nuevas_transiciones
        resultado.creacionDelta()
        resultado.q0 = f"{{{afd1.q0},{afd2.q0}}}"
        resultado.estadosInaccesibles = resultado.hallarEstadosInaccesibles()
        while resultado.estadosInaccesibles:
            x = resultado.estadosInaccesibles[0]
            resultado.delta = np.delete(resultado.delta, resultado.Q.index(x), 0)
            resultado.Q.remove(x)
            resultado.estadosInaccesibles = resultado.hallarEstadosInaccesibles()
        return resultado

    def hallarProductoCartesianoDiferenciaSimetrica(self, afd1, afd2):
        estados = list()
        estados_finales = list()
        nuevas_transiciones = list()
        resultado = copy.deepcopy(afd1)
        for estado_a1 in afd1.Q:
            for estado_a2 in afd2.Q:
                for elemento in resultado.Sigma:
                    nuevas_transiciones.append(
                        f"{{{estado_a1},{estado_a2}}}:{elemento}>{{{afd1.delta[afd1.Q.index(estado_a1), afd1.Sigma.index(elemento)]},"
                        f"{afd2.delta[afd2.Q.index(estado_a2), afd2.Sigma.index(elemento)]}}}")
                    print(
                        f"\u03B4(({estado_a1},{estado_a2}),{elemento}) = (\u03B4\N{SUBSCRIPT ONE}({estado_a1},{elemento}),\u03B4\N{SUBSCRIPT TWO}({estado_a2},{elemento}))="
                        f"({afd1.delta[afd1.Q.index(estado_a1), afd1.Sigma.index(elemento)]},{afd2.delta[afd2.Q.index(estado_a2), afd2.Sigma.index(elemento)]}),")
                estados.append(f"{{{estado_a1},{estado_a2}}}")
                if (estado_a1 in afd1.F and estado_a2 not in afd2.F) or (
                        estado_a1 not in afd1.F and estado_a2 in afd2.F): estados_finales.append(
                    f"{{{estado_a1},{estado_a2}}}")
        resultado.Q = estados
        resultado.F = estados_finales
        resultado.delta = nuevas_transiciones
        resultado.creacionDelta()
        resultado.q0 = f"{{{afd1.q0},{afd2.q0}}}"
        resultado.estadosInaccesibles = resultado.hallarEstadosInaccesibles()
        while resultado.estadosInaccesibles:
            x = resultado.estadosInaccesibles[0]
            resultado.delta = np.delete(resultado.delta, resultado.Q.index(x), 0)
            resultado.Q.remove(x)
            resultado.estadosInaccesibles = resultado.hallarEstadosInaccesibles()
        return resultado

    def hallarProductoCartesiano(self, afd1, afd2, operacion):
        if operacion == "interseccion":
            return self.hallarProductoCartesianoY(afd1, afd2)
        elif operacion == "union":
            return self.hallarProductoCartesianoO(afd1, afd2)
        elif operacion == "diferencia":
            return self.hallarProductoCartesianoDiferencia(afd1, afd2)
        elif operacion == "diferencia_simetrica":
            return self.hallarProductoCartesianoDiferenciaSimetrica(afd1, afd2)
        else:
            raise Exception(f"No existe el producto cartesiano de tipo {operacion}")

    def eliminarEstadosInaccesibles(self):
        self.estadosInaccesibles = self.hallarEstadosInaccesibles()
        while self.estadosInaccesibles:
            x = self.estadosInaccesibles[0]
            self.delta = np.delete(self.delta, self.Q.index(x), 0)
            self.Q.remove(x)
            self.estadosInaccesibles = self.hallarEstadosInaccesibles()

    def eliminarEstadosLimbo(self):
        self.estadosLimbo = self.hallarEstadosLimbo()
        while self.estadosLimbo:
            x = self.estadosLimbo[0]
            self.delta = np.delete(self.delta, self.Q.index(x), 0)
            self.Q.remove(x)
            self.estadosLimbo = self.hallarEstadosLimbo()

    def __str__(self):
        self.mostrarGrafoSimplificado("afd")
        return self.toString()
