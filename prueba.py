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
        #AFD que acepta cadenas que inician en 00
        a1=self.cargarAutomata(AFD,nombreArchivo="./automatas/afd.dfa")
        a1.mostrarGrafoSimplificado("./salidas/AFD/a1")
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
        a1.procesarListaCadenas(["01010101","10101010","11101010","00101010"],"./salidas/AFD/afda1_salida",imprimirPantalla=True)
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
        a2.procesarListaCadenas(["01010101","10101010","11101010","00001010"],"./salidas/AFD/a2_salida",imprimirPantalla=True)

    def probarAFN(self):
        an1=self.cargarAutomata(AFN,nombreArchivo="./automatas/afn.nfa")
        #Descomentar la siguiente línea para ver el grafo en imagen, GRAPHVIZ debe estar instalado.
        an1.mostrarGrafoSimplificado()
        print(an1.procesarCadena("aaaa"))
        print(an1.procesarCadena("aabb"))
        print(an1.procesarCadena("bbbb"))
        an1.procesarCadenaConDetalles("aaaa")
        an1.procesarCadenaConDetalles("aabb")
        an1.procesarCadenaConDetalles("bbbb")
        an1.computarTodosLosProcesamientos("aaaa","./salidas/AFN/afn1aaaa")
        an1.computarTodosLosProcesamientos("aabb","./salidas/AFN/afn1aabb")
        an1.computarTodosLosProcesamientos("bbbb","./salidas/AFN/afn1bbbb")

    def probarAFNLambda(self):
        an1 = self.cargarAutomata(AFNLambda, nombreArchivo="./automatas/afnl.nfe")
        #an1.mostrarGrafoSimplificado() #Mostrar grafo
        print(f"lambda clausura q0 -> {an1.lambdaClausura('q0')}")
        print(f"lambda clausura q0,q1,q2 -> {an1.calcularLambdaClausura('q0','q1','q2')}")
        print(an1.procesarCadena('bb')) #Procesar cadena a
        print(an1.procesarCadena('abb')) #Procesar cadena aba
        an1.procesarCadenaConDetalles('bb') #Procesar cadena bb con un único procesamiento de aceptación
        an1.procesarCadenaConDetalles('abb') #Procesar cadena aba con un único procesamiento de aceptación
        an1.computarTodosLosProcesamientos('abb',"./salidas/AFNL/an1abb")
        an1.procesarListaCadenas(['aba','bb'],"./salidas/AFNL/an1Lista",True)

    def probarAFNtoAFD(self):
        #AFN que acepta cadenas terminadas en b sobre {a,b}
        an1 = self.cargarAutomata(AFN, nombreArchivo="./automatas/afn1.nfa")
        an1.mostrarGrafoSimplificado("./salidas/AFNtoAFD/an1")
        af1 = an1.AFNtoAFD(an1)
        af1.mostrarGrafoSimplificado("./salidas/AFNtoAFD/af1")
        an1.procesarCadenaConDetalles("aaa")
        an1.procesarCadenaConDetalles("aaab")
        af1.procesarCadenaConDetalles("aaa")
        af1.procesarCadenaConDetalles("aaab")

    def probarAFNLambdaToAFN(self):
        # AFNLambda que acepta cadenas de la forma ab*a sobre {a,b}
        anl1 = self.cargarAutomata(AFNLambda, nombreArchivo="./automatas/afnl1.nfe")
        anl1.mostrarGrafoSimplificado("./salidas/AFNLambdaToAFN/anl1")
        an1 = anl1.AFN_LambdaToAFN(anl1)
        an1.mostrarGrafoSimplificado("./salidas/AFNLambdaToAFN/an1")
        anl1.procesarCadenaConDetalles("a")
        anl1.procesarCadenaConDetalles("aaa")
        anl1.procesarCadenaConDetalles("abbba")
        an1.procesarCadenaConDetalles("a")
        an1.procesarCadenaConDetalles("aaa")
        an1.procesarCadenaConDetalles("abbba")
        anl1.exportar("./salidas/AFNLambdaToAFN/anl1export")
        an1.exportar("./salidas/AFNLambdaToAFN/an1export")

    def probarAFNLambdaToAFD(self):
        # AFNLambda que acepta cadenas de la forma ab*a sobre {a,b}
        anl1 = self.cargarAutomata(AFNLambda, nombreArchivo="./automatas/afnl1.nfe")
        anl1.mostrarGrafoSimplificado("./salidas/AFNLambdaToAFD/anl1")
        an1 = anl1.AFN_LambdaToAFN(anl1)
        an1.mostrarGrafoSimplificado("./salidas/AFNLambdaToAFD/an1")
        ad1 = an1.AFNtoAFD(an1)
        ad1.mostrarGrafoSimplificado("./salidas/AFNLambdaToAFD/ad1")
        #NFE
        anl1.procesarCadenaConDetalles("a")
        anl1.procesarCadenaConDetalles("aaa")
        anl1.procesarCadenaConDetalles("abbba")
        #NFA
        an1.procesarCadenaConDetalles("a")
        an1.procesarCadenaConDetalles("aaa")
        an1.procesarCadenaConDetalles("abbba")
        #DFA
        ad1.procesarCadenaConDetalles("a")
        ad1.procesarCadenaConDetalles("aaa")
        ad1.procesarCadenaConDetalles("abbba")

        anl1.exportar("./salidas/AFNLambdaToAFD/anl1export")
        an1.exportar("./salidas/AFNLambdaToAFD/an1export")
        ad1.exportar("./salidas/AFNLambdaToAFD/ad1export")

    def probarComplemento(self):
        ad1 = self.cargarAutomata(AFD, nombreArchivo="./automatas/afd.dfa")
        ad1.mostrarGrafoSimplificado("./salidas/complemento/afdOriginal")
        ad2 = ad1.hallarComplemento(ad1)
        ad2.mostrarGrafoSimplificado("./salidas/complemento/afdComplemento")

    def probarProductoCartesiano(self):
        ad1 = self.cargarAutomata(AFD, nombreArchivo="./automatas/afd.dfa") #Lenguajes que empiezan por 00 sobre {0,1}
        ad2 = self.cargarAutomata(AFD, nombreArchivo="./automatas/afd1.dfa") #Lenguajes cuyo tercer elemento es un 0
        y = ad1.hallarProductoCartesiano(ad1,ad2,"interseccion") #Empiezan por 000
        y.mostrarGrafoSimplificado("./salidas/productoCartesiano/y")
        print(y.toString())
        o = ad1.hallarProductoCartesiano(ad1, ad2, "union")  #Empiezan por 00 o su tercer elemento es 0
        o.mostrarGrafoSimplificado("./salidas/productoCartesiano/o")
        print(o.toString())
        dif = ad1.hallarProductoCartesiano(ad1, ad2, "diferencia")  #Empiezan por 00 y su tercer elemento no es 0
        dif.mostrarGrafoSimplificado("./salidas/productoCartesiano/dif")
        print(dif.toString())
        sim = ad1.hallarProductoCartesiano(ad1, ad2, "diferencia_simetrica")  #No se cumplen ambas al tiempo
        sim.mostrarGrafoSimplificado("./salidas/productoCartesiano/sim")
        print(sim.toString())

    def probarSimplificacion(self):
        ad1 = self.cargarAutomata(AFD, nombreArchivo="./automatas/afd.dfa")  # Lenguajes que empiezan por 00 sobre {0,1}
        ad2 = self.cargarAutomata(AFD, nombreArchivo="./automatas/afd1.dfa")  # Lenguajes cuyo tercer elemento es un 0
        y = ad1.hallarProductoCartesiano(ad1, ad2, "interseccion")  # Empiezan por 000
        y_simplificado = y.simplificarAFD(y)
        y_simplificado.mostrarGrafoSimplificado("./salidas/simplificacion/y_simplificado")





    def main(self):
        #self.probarAFD()
        #self.probarAFN()
        #self.probarAFNLambda()
        #self.probarAFNtoAFD()
        #self.probarAFNLambdaToAFN()
        #self.probarAFNLambdaToAFD()
        #self.probarComplemento()
        #self.probarProductoCartesiano()
        #self.probarSimplificacion()



prueba = clasePrueba()
prueba.main()


