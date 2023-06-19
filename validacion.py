from afd import AFD
from afn import AFN
from afnl import AFNLambda
import random

class validacion:
    def cargarAutomata(self, tipo, **kwargs):
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

    def validarAFNtoAFD(self, *listaDeAFNs):
        """Por listaDeAFNs se entiende una lista de paths correspondientes a archivos de AFNs"""
        ans = [self.cargarAutomata(AFN, nombreArchivo=path) for path in listaDeAFNs]
        for i in range(len(ans)):
            print(f"Automata {i}:")
            automata = ans[i]
            afd_equivalente = automata.AFNtoAFD(automata)
            iguales = int()
            diferentes = int()
            cadenas_controversiales=list()
            for _ in range(5000):
                cadena = "".join([str(random.choice(ans[i].Sigma)) for _ in range(random.randint(50,100))])
                resultado_afn = automata.procesarCadena(cadena)
                resultado_afd = afd_equivalente.procesarCadena(cadena)
                if resultado_afn == resultado_afd: iguales+=1
                else:
                    cadenas_controversiales.append(cadena)
                    diferentes+=1
            print(f"Iguales: {iguales}. Diferentes: {diferentes}")
            print(f"cadenas que causaron diferencia: {cadenas_controversiales} \n")

    def validarAFNLambdaToAFN(self, *listaDeAFNLambdas):
        anls = [self.cargarAutomata(AFNLambda, nombreArchivo=path) for path in listaDeAFNLambdas]
        for i in range(len(anls)):
            print(f"Automata {i}:")
            automata = anls[i]
            afn_equivalente = automata.AFN_LambdaToAFN(automata)
            iguales = int()
            diferentes = int()
            cadenas_controversiales = list()
            for _ in range(5000):
                cadena = "".join([str(random.choice(afn_equivalente.Sigma)) for _ in range(random.randint(50, 100))])
                resultado_afnl = automata.procesarCadena(cadena)
                resultado_afn = afn_equivalente.procesarCadena(cadena)
                if resultado_afnl == resultado_afn:
                    iguales += 1
                else:
                    cadenas_controversiales.append(cadena)
                    diferentes += 1
            print(f"Iguales: {iguales}. Diferentes: {diferentes}")
            print(f"cadenas que causaron diferencia: {cadenas_controversiales} \n")



    def main(self):
        #self.validarAFNtoAFD("./automatas/afn.nfa", "./automatas/afn1.nfa")
        self.validarAFNLambdaToAFN("./automatas/afnl.nfe", "./automatas/afnl1.nfe")


valida = validacion()
valida.main()
