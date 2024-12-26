from qiskit import QuantumCircuit, transpile, QuantumRegister, ClassicalRegister
from qiskit_aer import AerSimulator
from  qiskit.quantum_info import Kraus
from qiskit.quantum_info import Statevector
from qiskit.quantum_info import Operator


iop=Operator.from_label('I')
xop=Operator.from_label('X')
yop=Operator.from_label('Y')
zop=Operator.from_label('Z')
hop=Operator.from_label('H')
sop=Operator.from_label('S')
proj0op=Operator.from_label('0')
proj1op=Operator.from_label('1')


#Use lambda function to define the K gate
Kop = lambda p: Operator([[1,0],[0,p**0.5]])
iexpand = lambda n: Operator.from_label('I'*n)


#Act gate on a qubit
single_op_on_qubit = lambda opstr,ntotal,qubitindex: Operator.from_label('I'*(qubitindex)) + opstr + Operator.from_label('I'*(ntotal-qubitindex))


proj_on_qubit = lambda projop,ntotal,qubitindex: (
                projop if ntotal==1 else
                projop ^ iexpand(ntotal-1) if qubitindex == 0 else
                iexpand(qubitindex) ^ projop if qubitindex == ntotal - 1 else
                iexpand(qubitindex) ^ projop ^ iexpand(ntotal-qubitindex-1)
)


Kop_onqubit = lambda p, ntotal, qubitindex: (
    Kop(p) if ntotal==1 else
    Kop(p) ^ iexpand(ntotal - qubitindex - 1) if qubitindex == 0 else
    iexpand(qubitindex) ^ Kop(p) if qubitindex == ntotal - 1 else
    iexpand(qubitindex) ^ Kop(p) ^ iexpand(ntotal - qubitindex - 1)
)


def cnot_on_qubits(ntotal,qindex1,qindex2):
    circ=QuantumCircuit(ntotal)
    circ.cx(ntotal-1-qindex1,ntotal-1-qindex2)
    return Operator.from_circuit(circ)




#Qiskit doesn't not support non-unitary circuit simulation
#Noisy Quantum Circuit which simulate a quantum circuit with non-unitary noise
#by naive state vector simulation
class noisyQcircuit:

    def __init__(self,n,t,p):
        self._nqubit=0
        self._t=t
        self._p=p
        self._circuit_description=[]#List of gates for the circuit
        self._initstate=Statevector.from_label('0'*n)
        self._state=None

    def add_x(self, qindex):
        self._circuit_description.append(('x',qindex))

    def add_y(self, qindex):
        self._circuit_description.append(('y',qindex))

    def add_z(self, qindex):
        self._circuit_description.append(('z',qindex))

    def add_s(self, qindex):
        self._circuit_description.append(('s',qindex))

    def add_cnot(self, qindex1, qindex2):
        self._circuit_description.append(('cnot',(qindex1,qindex2)))


    def add_hadamard(self, qindex):
        self._circuit_description.append(('h',qindex))


    #Add a projection gate to 0
    def add_proj0(self, qindex):
        self._circuit_description.append(('proj0',qindex))

    #Add a projection gate to 1
    def add_proj1(self, qindex):
        self._circuit_description.append(('proj1',qindex))

    #Inject a noise with Kgadget method
    def inject_noise(self, qindex):
        self._circuit_description.append(('K',qindex))
  

    #Run the circuit and return a final state vector
    def run(self):
        state=self._initstate 
        for (gate,qindex) in self._circuit_description:
            if gate=='x':
                state=state.evolve(single_op_on_qubit('X',self._nqubit,qindex))
            elif gate=='y':
                state=state.evolve(single_op_on_qubit('Y',self._nqubit,qindex))
            elif gate=='z':
                state=state.evolve(single_op_on_qubit('Z',self._nqubit,qindex))
            elif gate=='s':
                state=state.evolve(single_op_on_qubit('S',self._nqubit,qindex))
            elif gate=='h':
                state=state.evolve(single_op_on_qubit('H',self._nqubit,qindex))
            elif gate=='proj0':
                state=state.evolve(proj_on_qubit(proj0op,self._nqubit,qindex))
            elif gate=='proj1':
                state=state.evolve(proj_on_qubit(proj1op,self._nqubit,qindex))
            elif gate=='cnot':
                state=state.evolve(cnot_on_qubits(self._nqubit,qindex[0],qindex[1]))
            elif gate=='K':
                state=state.evolve(Kop_onqubit(self._p,self._nqubit,qindex))
        self._state=state
        return state





class KGadgetSimulator:


    def __init__(self,n,t,p):
        #Define Paramters
        self._n = n#Number of qubits
        self._t = t#Number of different noise
        self._p = p#Parameter of amplititude damping noise

        self._exact_noiseless_simulator = AerSimulator(method="stabilizer")#Qiskit stabilizer simulator for exact simulation without any noise
        self._exact_noiseless_circuit=QuantumCircuit(n,n)


        self._exact_noise_simulator = noisyQcircuit(n,t,p)#Qiskit stabilizer simulator for exact simulation with  noise


        self._kgadget_circ=None#Circuit for Kgadget simulation
        self._kgadget_simulator= AerSimulator(method="stabilizer")#Qiskit stabilizer simulator for Kgadget simulation


        self._circuit_description=[]#List of gates for the circuit
        self._Kraus_gateindex=[]#List of Kraus index for each noise gate

        self._current_gateindex=0#Index of the current gate


    def add_x(self, qindex):
        self._exact_noiseless_circuit.x(qindex)
        self._exact_noise_simulator.add_x(qindex)
        self._circuit_description.append(('x',qindex))
        self._current_gateindex+=1

    def add_y(self, qindex):
        self._exact_noiseless_circuit.y(qindex)
        self._exact_noise_simulator.add_y(qindex)
        self._circuit_description.append(('y',qindex))
        self._current_gateindex+=1

    def add_z(self, qindex):
        self._exact_noiseless_circuit.z(qindex)
        self._exact_noise_simulator.add_z(qindex)
        self._circuit_description.append(('z',qindex))
        self._current_gateindex+=1


    def add_s(self, qindex):
        self._exact_noiseless_circuit.s(qindex)
        self._exact_noise_simulator.add_s(qindex)
        self._circuit_description.append(('s',qindex))
        self._current_gateindex+=1

    def add_cnot(self, qindex1, qindex2):
        self._exact_noiseless_circuit.cx(qindex1,qindex2)
        self._exact_noise_simulator.add_cnot(qindex1,qindex2)
        self._circuit_description.append(('cnot',(qindex1,qindex2)))
        self._current_gateindex+=1


    def add_hadamard(self, qindex):
        self._exact_noiseless_circuit.h(qindex)
        self._exact_noise_simulator.add_hadamard(qindex)
        self._circuit_description.append(('h',qindex))
        self._current_gateindex+=1


    #Inject a noise with Kgadget method
    def inject_noise(self, qindex):
        self._circuit_description.append(('K',qindex))
        self._exact_noise_simulator.inject_noise(qindex)
        self._Kraus_gateindex.append(self._current_gateindex)
        self._current_gateindex+=1



    #Initialize the exact state vector of kgadget method
    def exact_kstates(self):
        pass



    def compile_kgadget_circuit(self):
        dataqubit=QuantumRegister(self._n,name='data')
        noisequbit=QuantumRegister(self._t,name='noise')
        self._kgadget_circ= QuantumCircuit(dataqubit,noisequbit)
        tmpKindex=0
        for (gate,qindex) in self._circuit_description:
            if gate=='x':
                self._kgadget_circ.x(qindex)
            elif gate=='y':
                self._kgadget_circ.y(qindex)
            elif gate=='z':
                self._kgadget_circ.z(qindex)
            elif gate=='s':
                self._kgadget_circ.s(qindex)
            elif gate=='h':
                self._kgadget_circ.h(qindex)
            elif gate=='cnot':
                self._kgadget_circ.cx(qindex[0],qindex[1])
            elif gate=='K':
                self._kgadget_circ.cx(qindex,self._n+tmpKindex)
                tmpKindex+=1


    def run_kadget_circuit(self,inputstr):
        simcirc=QuantumCircuit(self._n+self._t)
        for i in range(self._t):
            if inputstr[i]=='+':
                simcirc.h(self._n+i)
        simcirc.append(self._kgadget_circ,range(self._n))
        simcirc.save_statevector()

        stabilizer_simulator = AerSimulator(method="statevector")
        circ = transpile(circ, stabilizer_simulator)
        # Run and get statevector
        result = stabilizer_simulator.run(circ).result()
        statevector = result.get_statevector(simcirc)
        return statevector
        
        

    #Calculate the fidelity after run several simulations        
    def fidelity(self):
        pass


    #Run the exact circuit without any noise and return a final state vector
    def run_exact_noiseless(self):
        pass


    def run_exact_noisy(self):
        return self._exact_noise_simulator.run()