import matplotlib.pyplot as plt
from graphviz import Digraph,Source
import numpy as np
import copy
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
        if '-' in self.Sigma[0]:
            self.Sigma=[*self.Sigma,*self.rangoLenguaje(self.Sigma[0])]
            self.Sigma.pop(0)
            self.Sigma=sorted(self.Sigma)
        self.delta,self.Q=self.creacionDelta(self.delta,self.Sigma,self.Q)
        self.estadosLimbo=self.hallarEstadosLimbo(self.delta,self.Q,self.Sigma)
        self.estadosInaccesibles=self.hallarEstadosInaccesibles(self.delta,self.Q,self.q0)

    def rangoLenguaje(self,lenguaje):
        return [chr(caracter) for caracter in range(ord(lenguaje[0]),ord(lenguaje[-1])+1)]

    def manejarTransiciones(self,transiciones):
        return [re.split(':|>', x) for x in transiciones]

    def creacionDelta(self,transiciones,alfabeto,estados):
        matriz=np.asarray([['empt' for _ in range(len(alfabeto))] for _ in range(len(estados))],dtype=np.dtype('U100'))
        transiciones=self.manejarTransiciones(transiciones)
        for x in transiciones:
            if x[2] not in estados: raise Exception(f"La transición {x[0]}({x[1]}) lleva al estado {x[2]} que no existe entre los estados creados")
            matriz[estados.index(x[0]),alfabeto.index(x[1])]=x[2]
        matriz,estados=self.verificarCorregirCompletitudAFD(matriz,estados,alfabeto)
        return matriz,estados

    def verificarCorregirCompletitudAFD(self,matriz,estados,alfabeto):
        if 'empt' in matriz:
            estados.append('l')
            matriz[matriz=='empt']='l'
            matriz=np.vstack((matriz,['l' for _ in range(len(alfabeto))]))
        return matriz,estados

    def hallarEstadosLimbo(self,matriz,estados,alfabeto):
        estados_limbo=list()
        for i in range(len(matriz)):
            if np.array_equal(matriz[i],np.array(['l' for _ in range(len(alfabeto))])): estados_limbo.append(estados[i])
        return estados_limbo

    def hallarEstadosInaccesibles(self,matriz,estados,estado_inicial):
        estados_inaccesibles=list()
        for x in estados:
            if (x not in matriz or (len(np.where(matriz==x)[0].tolist())==1 and np.where(matriz==x)[0].tolist()==[estados.index(x)])) and x!=estado_inicial: estados_inaccesibles.append(x)
        return estados_inaccesibles

    def retornarTransiciones(self,matriz,estados,alfabeto):
        transiciones=list()
        for i in range(len(estados)):
            for j in range(len(alfabeto)):
                transiciones.append(f"{estados[i]}:{alfabeto[j]}>{matriz[i,j]}")
        return transiciones

    def imprimirAFDSimplificado(self):
        self.mostrarGrafoSimplificado()
        inaccesibles=self.hallarEstadosInaccesibles(self.delta,self.Q,self.q0)
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
        self.mostrarGrafo()
        return "#!dfa\n"+"#alphabet\n"+"\n".join(self.Sigma)+"\n#states\n"+"\n".join(self.Q)+"\n#initial\n"+f"{self.q0}\n"+ \
        "#accepting\n"+"\n".join(self.F)+"\n#transitions\n"+"\n".join(self.retornarTransiciones(self.delta,self.Q,self.Sigma))

    def mostrarGrafoSimplificado(self,archivo="grafo_s"):
        transiciones=self.manejarTransiciones(self.retornarTransiciones(self.delta,self.Q,self.Sigma))
        g_s=Digraph(strict=False)
        g_s.attr(rankdir="LR")
        g_s.node('i',label='',shape="point")
        g_s.edge('i',self.q0,arrowsize='0.5')
        for x in self.Q:
            if x in self.F:
                g_s.node(x,shape="doublecircle")
            elif x!='l' and x not in self.estadosInaccesibles:
                g_s.node(x,shape="circle")
        for x in transiciones:
            if 'l' not in x and x[0] not in self.estadosInaccesibles and x[2] not in self.estadosInaccesibles: g_s.edge(x[0],x[2],label=x[1],arrowsize='0.5')
        s=Source(g_s.source,filename=archivo,format="svg")
        s.view()

    def mostrarGrafo(self,archivo="grafo"):
        transiciones=self.manejarTransiciones(self.retornarTransiciones(self.delta,self.Q,self.Sigma))
        g=Digraph(strict=False)
        g.attr(rankdir="LR")
        g.node('i',label='',shape="point")
        g.edge('i',self.q0,arrowsize='0.5')
        for x in self.Q:
            if x in self.F: g.node(x,shape="doublecircle")
            else: g.node(x,shape="circle")
        for x in transiciones:
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

    def hallarProductoCartesianoY(self,afd1,afd2):
        estados=list()
        estados_finales=list()
        nuevas_transiciones=list()
        resultado=copy.deepcopy(afd1)
        for estado_a1 in afd1.Q:
            for estado_a2 in afd2.Q:
                for elemento in resultado.Sigma:
                    nuevas_transiciones.append(f"{{{estado_a1},{estado_a2}}}:{elemento}>{{{afd1.delta[afd1.Q.index(estado_a1),afd1.Sigma.index(elemento)]},"
                                               f"{afd2.delta[afd2.Q.index(estado_a2),afd2.Sigma.index(elemento)]}}}")
                    print(f"\u03B4(({estado_a1},{estado_a2}),{elemento}) = (\u03B4\N{SUBSCRIPT ONE}({estado_a1},{elemento}),\u03B4\N{SUBSCRIPT TWO}({estado_a2},{elemento}))="
                          f"({afd1.delta[afd1.Q.index(estado_a1),afd1.Sigma.index(elemento)]},{afd2.delta[afd2.Q.index(estado_a2),afd2.Sigma.index(elemento)]}),")
                estados.append(f"{{{estado_a1},{estado_a2}}}")
                if estado_a1 in afd1.F and estado_a2 in afd2.F: estados_finales.append(f"{{{estado_a1},{estado_a2}}}")
        resultado.Q=estados
        resultado.F=estados_finales
        resultado.delta=nuevas_transiciones
        resultado.delta,resultado.Q=resultado.creacionDelta(resultado.delta,resultado.Sigma,resultado.Q)
        resultado.q0=f"{{{afd1.q0},{afd2.q0}}}"
        resultado.estadosInaccesibles=resultado.hallarEstadosInaccesibles(resultado.delta,resultado.Q,resultado.q0)
        while resultado.estadosInaccesibles:
            x=resultado.estadosInaccesibles[0]
            resultado.delta=np.delete(resultado.delta,resultado.Q.index(x),0)
            resultado.Q.remove(x)
            resultado.estadosInaccesibles=resultado.hallarEstadosInaccesibles(resultado.delta,resultado.Q,resultado.q0)
        return resultado

    def __str__(self): return self.toString()

a=AFD('afd.dfa')
b=AFD('afd1.dfa')
c=a.hallarProductoCartesianoY(a,b)