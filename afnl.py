import matplotlib.pyplot as plt
from graphviz import Digraph, Source
from collections import deque
from afd import AFD
from afn import AFN
import numpy as np
import copy
import re


class AFNLambda(AFN):
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

    def ANF_LambdaToAFN(self,afnl):
        afn=copy.deepcopy(afnl)
        transiciones=afn.delta
        transiciones=np.delete(transiciones,-1,1)
        for estado in afn.Q:
            clausuras=afnl.calcularLambdaClausura(estado)
            if len(set(afn.F).intersection(set(np.array(clausuras).flatten())))>0: afn.F.append(estado)
            for i in range(len(afn.Sigma)-1):
                nueva_transicion = list()
                for clausura in clausuras:
                    if 'l' not in afn.delta[afn.Q.index(clausura),i]: nueva_transicion.append(afn.delta[afn.Q.index(clausura),i])
                if len(nueva_transicion)==0: transiciones[afn.Q.index(estado),i]='l'
                else:
                    transiciones[afn.Q.index(estado),i]=';'.join(set(np.array(afn.calcularLambdaClausura(*nueva_transicion)).flatten()))
        return AFN(alfabeto=afn.Sigma[:-1],estados=afn.Q,estadoInicial=afn.q0,estadosAceptacion=afn.F,transiciones=transiciones,deltaEnFormato=True)



    def __str__(self):
        return "#!nfe" + self.toString()[5:]



a = AFNLambda(nombreArchivo="afnl1.nfe",esLambda=True)
b=a.ANF_LambdaToAFN(a)
b.imprimirAFNSimplificado()