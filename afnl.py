import matplotlib.pyplot as plt
from graphviz import Digraph, Source
from collections import deque
from afd import AFD
from afn import AFN
import numpy as np
import itertools
import copy
import re


class AFNLambda(AFN):
    def __init__(self,nombreArchivo):
        super().__init__(nombreArchivo=nombreArchivo,esLambda=True)

    def lambdaClausura(self,estado):
        stack = [estado]
        clausura = {estado}
        while stack:
            current_state = stack.pop()
            if ';' in self.delta[self.Q.index(current_state),-1]:
                epsilon_transitions = self.delta[self.Q.index(current_state),-1].split(';')
                for transition in epsilon_transitions:
                    if transition not in clausura and 'l' not in transition:
                        clausura.add(transition)
                        stack.append(transition)
            else:
                if self.delta[self.Q.index(current_state),-1] not in clausura and 'l' not in self.delta[self.Q.index(current_state),-1]:
                    clausura.add(self.delta[self.Q.index(current_state),-1])
                    stack.append(self.delta[self.Q.index(current_state),-1])

        return sorted(clausura)

    def calcularLambdaClausura(self,*estados):
        clausuras=list()
        for estado in estados:
            clausuras.append(self.lambdaClausura(estado))
        if len(clausuras)>1: return clausuras
        return clausuras[0]

    def imprimirANFLSimplificado(self):
        self.imprimirAFNSimplificado()

    def exportar(self,archivo):
        contenido = "#!nfe" + self.toString()[5:]
        with open(f"{archivo}.nfe", 'w') as f:
            f.write(contenido)

    def AFN_LambdaToAFN(self,afnl,imprimir=True):
        afn=copy.deepcopy(afnl)
        transiciones=afn.delta
        transiciones=np.delete(transiciones,-1,1)
        for estado in afn.Q:
            clausuras=afnl.calcularLambdaClausura(estado)
            if imprimir: print(f"\u03BB[{estado}]={{{','.join(clausuras)}}}")
            if len(set(afn.F).intersection(set(np.array(clausuras).flatten())))>0: afn.F.append(estado)
            for i in range(len(afn.Sigma)-1):
                nueva_transicion = list()
                if imprimir: print(f"\u0394({estado},{afn.Sigma[i]}) = \u03BB[\u0394(\u03BB[{estado}],{afn.Sigma[i]})] = "
                      f"\u03BB[\u0394({{{clausuras}}},{afn.Sigma[i]})]",end=" = ")
                for clausura in clausuras:
                    if ';' in afn.delta[afn.Q.index(clausura),i]: nueva_transicion=[*nueva_transicion,*afn.delta[afn.Q.index(clausura),i].split(';')]
                    elif 'l' not in afn.delta[afn.Q.index(clausura),i]: nueva_transicion.append(afn.delta[afn.Q.index(clausura),i])
                if imprimir: print(f"\u03BB[{{{nueva_transicion}}}]",end=" = ")
                if len(nueva_transicion)==0: transiciones[afn.Q.index(estado),i]='l'
                elif len(nueva_transicion)==1: transiciones[afn.Q.index(estado),i]=';'.join(set(afn.calcularLambdaClausura(*nueva_transicion)))
                else: transiciones[afn.Q.index(estado),i]=';'.join(set(np.concatenate(afn.calcularLambdaClausura(*nueva_transicion))))
                if imprimir: print(f"{transiciones[afn.Q.index(estado),i]}")
        return AFN(alfabeto=afn.Sigma[:-1],estados=afn.Q,estadoInicial=afn.q0,estadosAceptacion=afn.F,transiciones=transiciones,deltaEnFormato=True)

    def AFN_LambdaToAFD(self,afnl):
        afn=self.ANF_LambdaToAFN(afnl)
        afd=afn.AFNtoAFD(afn)
        return afd

    def procesarCadena(self,cadena):
        alphabet = self.Sigma
        if '$' in alphabet: alphabet.remove('$')
        states = self.Q
        current_states = {self.q0}
        for char in cadena:
            next_states = set()
            for state in current_states:
                col_idx = alphabet.index(char)
                next_state_str = self.delta[states.index(state)][col_idx]
                next_states.update(next_state_str.split(';'))
            current_states = next_states
        while True:
            next_states = set()
            for state in current_states:
                eps_idx = len(alphabet)
                next_state_str = self.delta[states.index(state)][eps_idx]
                next_states.update(next_state_str.split(';'))
            if not next_states - current_states:
                break
            current_states.update(next_states)
        final_states = {x for x in self.F}
        if current_states.intersection(final_states):
            return True
        else:
            return False

    def procesarCadenaConDetalles(self, cadena):
        alphabet = copy.deepcopy(self.Sigma)
        states = copy.deepcopy(self.Q)
        initial_state = copy.deepcopy(self.q0)
        final_states = {x for x in self.F}
        if '$' in alphabet: alphabet.remove('$')
        stack = [(initial_state, cadena, [])]
        while stack:
            current_state, remaining_string, path = stack.pop()
            if not remaining_string:
                if current_state in final_states:
                    path.append(f"[{current_state},{remaining_string}]")
                    accepting_path = "->".join(path)
                    print(f"{accepting_path} -> Aceptación")
                    return True
                else:
                    continue
            char = remaining_string[0]
            col_idx = alphabet.index(char)
            next_state_str = self.delta[states.index(current_state)][col_idx]
            next_states = next_state_str.split(';')
            for state in next_states:
                stack.append((state, remaining_string[1:], path + [f"[{current_state},{remaining_string}]"]))
        print("No Aceptación")
        return False

    def computarTodosLosProcesamientos(self,cadena,nombreArchivo):
        a = copy.deepcopy(self)
        b = a.AFN_LambdaToAFN(a, False)
        return b.computarTodosLosProcesamientos(cadena,nombreArchivo)

    def procesarListaCadenas(self,listaCadenas,nombreArchivo,imprimirPantalla):
        a = copy.deepcopy(self)
        b = a.AFN_LambdaToAFN(a, False)
        return b.procesarListaCadenas(listaCadenas,nombreArchivo,imprimirPantalla)

    def procesarCadenaConversion(self,cadena):
        a = copy.deepcopy(self)
        b = a.AFN_LambdaToAFN(a, False)
        c = b.AFNtoAFD(b)
        return c.procesarCadena(cadena)

    def procesarCadenaConDetallesConversion(self,cadena):
        a = copy.deepcopy(self)
        b = a.AFN_LambdaToAFN(a, False)
        c = b.AFNtoAFD(b)
        return c.procesarCadenaConDetalles(cadena)

    def procesarListaCadenasConversion(self,listaCadenas,nombreArchivo,imprimirPantalla):
        a = copy.deepcopy(self)
        b = a.AFN_LambdaToAFN(a, False)
        c = b.AFNtoAFD(b)
        return c.procesarListaCadenas(listaCadenas,nombreArchivo,imprimirPantalla)

    def __str__(self):
        self.mostrarGrafoSimplificado("afnl")
        return "#!nfe" + self.toString()[5:]