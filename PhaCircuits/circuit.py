#circuit.py
#After import only the necessary parts of the librariesimport schemdraw
import schemdraw
import schemdraw.elements as elm
import numpy as np
from math import atan
import cmath

class Circuit:
    def __init__(self):
        self.elements =0 # Counts the number of elements
        self.elements_list =[] #List that saves the elements
        self.nodes = set()   #Set() makes the nodes stored in the set unique
        self.node_map={}     
        self.reference_node=0  #Default reference node
        self.currents ={}
        self.vsource_wire={}
        self.inches_punit()


    #Add Element
    def element(self, Element: str, Start: tuple, End: tuple, Value: float=0, Label: str = None):
        self.elements_list.append([Element, Start, End, Value, Label])
        #Add nodes to the set
        if Element not in 'Coupling':
            self.elements +=1
            self.nodes.add(Start)
            self.nodes.add(End)

    def inches_punit(self,value: float = 1.2):
        schemdraw.config(lw=1, font='serif',inches_per_unit=value)


    def show_elements(self):
        for i, element in enumerate(self.elements_list):
            print(f"Element {i+1}: {element}")


    def draw(self):
        with schemdraw.Drawing() as d:
            for iten in self.elements_list:
                if(iten[0]=='Resistor'):
                            d.add(elm.Resistor(scale=0.5).endpoints(iten[1], iten[2]).label( str(iten[3])+" Ω",loc='top', fontsize=8).label(iten[4],loc='bot', fontsize=8))
                elif(iten[0]=='Capacitor'):
                        capacitor_value = f"{int(iten[3].imag)}j" if iten[3].imag % 1 == 0 else f"{iten[3].imag}j"
                        C= d.add(elm.Capacitor(scale=0.7).endpoints(iten[1], iten[2]).label(capacitor_value + " Ω", loc='top', fontsize=8).label(iten[4], loc='bot', fontsize=8))
                elif(iten[0]=='Inductor'):
                        d.add(elm.Inductor(scale=0.7).endpoints(iten[1], iten[2]).label(str(iten[3])+" Ω",loc='top', fontsize=8).label(iten[4],loc='bot', fontsize=8))
                elif(iten[0] == 'Voltage Source'):
                        d.add(elm.SourceV(scale=0.7).endpoints(iten[1], iten[2]).label(str(iten[3])+"V",loc='top', fontsize=8).label(iten[4],loc='bot', fontsize=8)),
                elif(iten[0] == 'Current Source'):
                        d.add(elm.SourceI(scale=0.7).endpoints(iten[1], iten[2]).label(str(iten[3])+"A",loc='top', fontsize=8).label(iten[4],loc='bot', fontsize=8))
                elif(iten[0] == 'Wire'):
                        d.add(elm.Line(scale=0.5).endpoints(iten[1], iten[2]))

    def map_nodes(self):
        #Create a mapping to associeate each node with an index
        self.node_map = {}
        #Populate node_map dictionary
        # enumerate returns index and node 
        for idx, node in enumerate(self.nodes):
            self.node_map[node] = idx #The node is the key and idx is the value
        print("Mapeamento de nós:", self.node_map)


    def build_equations(self):
        self.map_nodes()
        n = len(self.nodes)
        G = np.zeros((n, n), dtype=complex) #Conductance Matrix
        I = np.zeros(n, dtype=complex) #Currents Matrix

        for elm in self.elements_list:
            Element, Start, End, Value, Label = elm

            if(Element == 'Coupling'):
                idxes_start = list()
                idxes_end = list()
                impedances = list()
                name1, name2, k = elm[1], elm[2], elm[3]

                for x in self.elements_list:
                    if x[-1] == name1 or x[-1] == name2:
                        Start, End = x[1], x[2]
                        idxes_start.append(self.node_map[Start])
                        idxes_end.append(self.node_map[End])
                        impedances.append(x[-2])

                coupling_impedance = k*np.sqrt(abs(impedances[0]*impedances[1]))*1j

                for i in range(2):
                    G[idxes_start[i], idxes_start[i]] -= 1 / impedances[i]
                    G[idxes_end[i], idxes_end[i]] -= 1 / impedances[i]
                    G[idxes_start[i], idxes_end[i]] += 1 / impedances[i]
                    G[idxes_end[i], idxes_start[i]] += 1 / impedances[i]

                    G[idxes_start[i], idxes_start[i]] += impedances[1]/(impedances[0]*impedances[1]*(1-k**2))
                    G[idxes_end[i], idxes_end[i]] += impedances[0]/(impedances[0]*impedances[1]*(1-k**2))
                    G[idxes_start[i], idxes_end[i]] -= coupling_impedance/(impedances[0]*impedances[1]*(1-k**2))
                    G[idxes_end[i], idxes_start[i]] -= coupling_impedance/(impedances[0]*impedances[1]*(1-k**2))
            else:
                idx_start=self.node_map[Start]
                idx_end=self.node_map[End]

                if(Element=='Current Source'):
                    I[idx_start] -= Value
                    I[idx_end] += Value

                #Voltage Source stamp from modified nodal analysis
                elif(Element == 'Voltage Source' ):
                    G = np.vstack([G, np.zeros(G.shape[1])])  #Add a new row
                    G = np.hstack([G, np.zeros((G.shape[0], 1))])  #Add a new column
                    G[-1, idx_start] = -1
                    G[-1, idx_end] = +1
                    G[idx_start,-1] =-1
                    G[idx_end,-1] =+1
                    I = np.append(I, Value)  
                    self.vsource_wire[Label] = I.shape[0]

                #Wire is a Voltage Source with valeu null
                elif(Element == 'Wire'):
                    G = np.vstack([G, np.zeros(G.shape[1])])  
                    G = np.hstack([G, np.zeros((G.shape[0], 1))])  
                    G[-1, idx_start] = -1
                    G[-1, idx_end] = +1
                    G[idx_start,-1] =-1
                    G[idx_end,-1] =+1
                    I = np.append(I, 0) 
                    
                #In AC analysis Resistor, Inductor and Capacitor is measured in ohms
                #So uses the resistor stamp
                elif(Element == 'Resistor' or Element == 'Inductor' or Element == 'Capacitor' ):
                    G[idx_start, idx_start] += 1 / Value
                    G[idx_end, idx_end] += 1 / Value
                    G[idx_start, idx_end] -= 1 / Value
                    G[idx_end, idx_start] -= 1 / Value

        return G, I

    def solve_circuit(self):
        G, I = self.build_equations()
        ref_node = self.reference_node

        # Reduce the matrix by removing the row and column from the reference node
        G_reduced = np.delete(np.delete(G, ref_node, axis=0), ref_node, axis=1)
        I_reduced = np.delete(I, ref_node)

        try:
            # Solve the system G * V = I
            V_reduced = np.linalg.solve(G_reduced, I_reduced)

            # Reconstruct full vector of voltage 
            self.V = np.insert(V_reduced, ref_node, 0)  # Define voltage in reference as zero
            #print("Tensões nos nós:", self.V)
        except np.linalg.LinAlgError as e:
            print("Error resolving system:", e)
            self.V = None
        return self.V

    def set_reference(self,reference: int):
        self.reference_node = reference

    
    def show_conductance_matrix(self):
        G, I = self.build_equations()
        row,colunm = G.shape
        for i in range (row):
          for j in range (colunm):
            print("{:.3f}".format(G[i,j]), "|",end=" ")
          print(" ")

    def calculate_currents(self):
        for elm in self.elements_list:
            Element, Start, End, Value, Label = elm
            idx_start=self.node_map[Start]
            idx_end=self.node_map[End]
            if Element =='Resistor' or Element =="Capacitor"  or  Element =="Inductor":
                self.currents[Label] = (self.V[idx_start]-self.V[idx_end])/Value
            elif Element == 'Voltage Source':
                indice = self.vsource_wire[Label]
                self.currents[Label] =  self.V[indice-1]
        

    def get_current(self, Label: str):
        self.calculate_currents()
        #print("Corrente ",self.currents[Label])
        return(self.currents[Label])

    def current_complex(self,value):
        modulo = abs(value)
        angulo = cmath.phase(value) * 180/np.pi
        return(modulo,angulo)
    
    def draw_with_currents(self,Label: list = [None]):

            with schemdraw.Drawing() as d:
                for iten in self.elements_list:
    
                    if(iten[0]=='Resistor'):
                        modulo,angulo=self.current_complex(self.get_current(iten[4]))
                        R=d.add(elm.Resistor(scale = 0.5).endpoints(iten[1], iten[2]).label( str(iten[3])+" Ω",loc='top', fontsize=8))
                        for name in Label:
                            if(iten[4]==name or name==None):
                                d.add(elm.CurrentLabelInline(direction='in',headlength=0.15,headwidth=0.15,ofst=0.5).at(R).label("{:.2f}".format(modulo)+"∠"+ "{:.2f}".format(angulo)+" A",fontsize=8  ) )
                                # d.add(elm.CurrentLabel(length=0.7,top=False).at(R).label("{:.2f}".format(modulo)+"∠"+ "{:.2f}".format(angulo)+" A",fontsize=8  ) )
                    elif(iten[0]=='Capacitor'):
                        modulo,angulo=self.current_complex(self.get_current(iten[4]))
                        capacitor_value = f"{int(iten[3].imag)}j" if iten[3].imag % 1 == 0 else f"{iten[3].imag}j"
                        C= d.add(elm.Capacitor(scale = 0.7).endpoints(iten[1], iten[2]).label(capacitor_value + " Ω", loc='top', fontsize=8))
                        for name in Label:
                            if(iten[4]==name or name==None):
                                d.add(elm.CurrentLabelInline(direction='in',headlength=0.15,headwidth=0.15,ofst=0.5).at(C).label("{:.2f}".format(modulo)+"∠"+ "{:.2f}".format(angulo)+" A",fontsize=8  ) )
                                # d.add(elm.CurrentLabel(length=0.7,top=False).at(C).label("{:.2f}".format(modulo)+"∠"+ "{:.2f}".format(angulo)+" A",fontsize=8  ) )
                    elif(iten[0]=='Inductor'):
                        modulo,angulo=self.current_complex(self.get_current(iten[4]))
                        L = d.add(elm.Inductor(scale = 0.7).endpoints(iten[1], iten[2]).label(str(iten[3])+" Ω",loc='top', fontsize=8))
                        for name in Label:
                            if(iten[4]==name or name==None):
                                d.add(elm.CurrentLabelInline(direction='in',headlength=0.15,headwidth=0.15,ofst=0.5).at(L).label("{:.2f}".format(modulo)+"∠"+ "{:.2f}".format(angulo)+" A",fontsize=8  ) )
                                # d.add(elm.CurrentLabel(length=0.7,top=False).at(L).label("{:.2f}".format(modulo)+"∠"+ "{:.2f}".format(angulo)+" A",fontsize=8  ) )
                    elif(iten[0] == 'Voltage Source'):
                        modulo,angulo=self.current_complex(self.get_current(iten[4]))
                        V=d.add(elm.SourceV(scale = 0.7).endpoints(iten[1], iten[2]).label(str(iten[3])+"V",loc='top', fontsize=8))
                        for name in Label:
                            if(iten[4]==name or name==None):
                                d.add(elm.CurrentLabelInline(direction='out',headlength=0.15,headwidth=0.15,ofst=0.5).at(V).label("{:.2f}".format(modulo)+"∠"+ "{:.2f}".format(angulo)+" A",fontsize=8  ) )
                                # d.add(elm.CurrentLabel(length=0.7,top=False,reverse=True).at(V).label("{:.2f}".format(modulo)+"∠"+ "{:.2f}".format(angulo)+" A",fontsize=8  ) )
                    elif(iten[0] == 'Current Source'):
                        d.add(elm.SourceI(scale = 0.7).endpoints(iten[1], iten[2]).label(str(iten[3])+"A",loc='top', fontsize=8))
                    elif(iten[0] == 'Wire'):
                        d.add(elm.Line(scale = 0.5).endpoints(iten[1], iten[2]))
                        #d.add(elm.Ground).endpoints(iten[1],iten[2])