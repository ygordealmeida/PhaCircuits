# Summary

- [Phasorial Circuit Simulator](#phasorial-circuit-simulator)
- [Requirements](#requirements)
- [How To Use](#how-to-use)



<h1 align="center">Phasorial Circuit Simulator</h1>

  AC circuits are solved in the phasor domain, and so when we are studying or solving questions, whether from a book or in the classroom, it is common to leave aside the temporal domain and use circuits that are already given in their phasor format, without having to transform the components.

However, in most simulators, the answers are presented in temporal form, so if you want to compare them with your answers obtained in the phasor format, it is necessary to make a relationship between the components and the quantities.

The objective of this repository is to provide a simulator that is viable for simulating circuits directly in phasor form.

<h1 align="center">Simulador de Circuitos Fasoriais</h1>

  Circuitos CA são resolvidos no domínio fasorial, e por isso quando estamos estudando ou resolvendo questões seja de um livro ou de sala de aula é comum deixarmos de lado o domínio temporal e utilizarmos circuitos que ja estão dados em seu formato fasorial, sem que seja necessário fazer a transformação dos componentes.

  No entanto, na maioria dos simuladores, as respostas são dispostas na forma temporal, então se você deseja comparar com suas respostas obtidas no formato fasorial, é necessário fazer uma relação entre os componentes e as grandezas. 

  O objetivo desse repositório é dispor um simulador que seja viável simular circuitos diretamente na forma fasorial.

 <div align="center">
  <img src="https://github.com/user-attachments/assets/bfb27055-f661-476e-89e5-5719d1f4a3a9" alt="Circuit Image" width="500"/>
</div>

  ## Requirements

-  python>=3.10
-  numpy>=1.24.1
-  schemdraw>=0.19

  ## How to Install

  PhasorialCircuit can be installed from pip using

   ```
    pip install PhasorialCircuit
   ```

  ## How To Use

  First you need to have the PhaCircuits folder in your code directory
  
  Then import Circuit from PhaCircuits and instantiate a circuit
  
   ```
    from PhaCircuits import Circuit

    circuit1 = Circuit()
   ```
  Then you can place the components of youy circuit using:
   ```
  circuit1.element("Name of element", (X_start, Y_start), (X_end, Y_end), Value, "Component Label")
   ```
  When your circuit is finished you can see the diagram with:
   ```
    circuit1.draw()
   ```
  And see the currents in the components with:
  ```
    circuit1.see_currents("Component1 Label", "Component2 Label")
  ```
  if you want see all currents use Label = None
  
  For more see examples and documentation.
