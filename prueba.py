from afd import AFD
from afn import AFN
from afnl import AFNLambda
import numpy as np

class clasePrueba:

    def cargarAutomata(self,tipo,**kwargs):
        if "nombreArchivo" in kwargs:
            return tipo(nombreArchivo=kwargs.get("nombreArchivo"))
        else:
            alfabeto=kwargs.get("alfabeto")
            estados=kwargs.get("estados")
            estadoInicial = kwargs.get("estadoInicial")
            estadosAceptacion = kwargs.get("estadosAceptacion")
            transiciones = kwargs.get("transiciones")
            if isinstance(transiciones,np.ndarray): return tipo(alfabeto = alfabeto, estados = estados, estadoInicial = estadoInicial,
            estadosAceptacion = estadosAceptacion, transiciones = transiciones, deltaEnFormato = True)
            else: return tipo(alfabeto = alfabeto, estados = estados, estadoInicial = estadoInicial,
            estadosAceptacion = estadosAceptacion, transiciones = transiciones)

    def probarAFD(self):
        """
        Se crean AFD a partir de archivos almacenados en el directorio /automatas, el primero (a1) corresponde al AFD
        que acepta cadenas que empiezan en 00 sobre el alfabeto {0,1}, el segundo (a2) a cadenas cuyo tercer elemento
        es un cero.

        La función print(afd) muestra el GRAFO obtenido en un archivo de imagen,
        dicho archivo se abre de manera automática (SE REQUIERE QUE GRAPHVIZ ESTÉ INSTALADO EN LA MÁQUINA Y EN pip),
        un ejemplo de la salida para el primer afd creado aquí, se encuentra en el directorio /automatas/grafo.svg
        afd.imprimirAFDSimplificado hace lo mismo pero sin estados limbo o inaccesibles.
        ¡ESTO FUNCIONA PARA AFD,AFN y AFNL. Al igual que para los resultados de cualquier operación que retorne un automata,
        como productos cartesianos, conversiones, etc.!
        :return:
        """

        #AFD que acepta cadenas que inician en 00
        a1=self.cargarAutomata(AFD,nombreArchivo="./automatas/afd.dfa")
        #print(a1)
        print("Procesamiento a1: ")
        #Procesamiento cadenas a1
        print(a1.procesarCadena("01010101"))
        print(a1.procesarCadena("10101010"))
        print(a1.procesarCadena("11101010"))
        print(a1.procesarCadena("00101010"))
        #Procesamiento cadenas con detalle a1
        a1.procesarCadenaConDetalles("01010101")
        a1.procesarCadenaConDetalles("10101010")
        a1.procesarCadenaConDetalles("11101010")
        a1.procesarCadenaConDetalles("00101010")
        #Procesar lista de cadenas y guardar en el archivo a1_salida en el directorio salidas
        a1.procesarListaCadenas(["01010101","10101010","11101010","00101010"],"./salidas/a1_salida",imprimirPantalla=True)
        print("---")

        # AFD que acepta cadenas que tienen como tercer elemento un 0
        a2=self.cargarAutomata(AFD, nombreArchivo="./automatas/afd1.dfa")
        print("Procesamiento a2")
        # Procesamiento cadenas a2
        print(a2.procesarCadena("01010101"))
        print(a2.procesarCadena("10101010"))
        print(a2.procesarCadena("11101010"))
        print(a2.procesarCadena("00001010"))
        # Procesamiento cadenas con detalle a2
        a2.procesarCadenaConDetalles("01010101")
        a2.procesarCadenaConDetalles("10101010")
        a2.procesarCadenaConDetalles("11101010")
        a2.procesarCadenaConDetalles("00001010")
        # Procesar lista de cadenas y guardar en el archivo a1_salida en el directorio salidas
        a2.procesarListaCadenas(["01010101","10101010","11101010","00001010"],"./salidas/a2_salida",imprimirPantalla=True)

    def probarAFN(self):
        an1=self.cargarAutomata(AFN,nombreArchivo="./automatas/afn.nfa")
        #Descomentar la siguiente línea para ver el grafo en imagen, GRAPHVIZ debe estar instalado.
        #an1.imprimirAFNSimplificado()
        print(an1.procesarCadena("aaaa"))
        print(an1.procesarCadena("aabb"))
        print(an1.procesarCadena("bbbb"))
        an1.procesarCadenaConDetalles("aaaa")
        an1.procesarCadenaConDetalles("aabb")
        an1.procesarCadenaConDetalles("bbbb")
        an1.computarTodosLosProcesamientos("aaaa","./salidas/afn1aaaa")
        an1.computarTodosLosProcesamientos("aabb","./salidas/afn1aabb")
        an1.computarTodosLosProcesamientos("bbbb","./salidas/afn1bbbb")





a=clasePrueba()
#a.probarAFD()
a.probarAFN()