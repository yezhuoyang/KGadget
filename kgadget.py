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
setop=Operator.from_label('S')
proj0op=Operator.from_label('0')
proj1op=Operator.from_label('1')


#Use lambda function to define the K gate
Kop = lambda p: Operator([[1,0],[0,p**0.5]])
iexpand = lambda n: Operator.from_label('I'*n)
Kop_onqubit = lambda p,ntotal,qubitindex: iexpand(qubitindex-1) ^ Kop(p) ^ iexpand(ntotal-qubitindex-1) 





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
        pass

    #Add a projection gate to 1
    def add_proj1(self, qindex):
        pass

    #Inject a noise with Kgadget method
    def inject_noise(self, qindex):
        pass
  
    #Run the circuit and return a final state vector
    def run(self):
        pass



class KGadgetSimulator:


    def __init__(self,n,t,p):
        #Define Paramters
        self._n = n#Number of qubits
        self._t = t#Number of different noise
        self._p = p#Parameter of amplititude damping noise

        self._exact_noiseless_simulator = AerSimulator(method="stabilizer")#Qiskit stabilizer simulator for exact simulation without any noise
        self._exact_noiseless_circuit=QuantumCircuit(n,n)

        self._exact_noise_simulator = AerSimulator()#Qiskit stabilizer simulator for exact simulation with  noise
        self._exact_noisy_circuit=QuantumCircuit(n,n)


        self._kgadget_circ=None#Circuit for Kgadget simulation
        self._kgadget_simulator= AerSimulator(method="stabilizer")#Qiskit stabilizer simulator for Kgadget simulation

        self._circuit_description=[]#List of gates for the circuit
        self._Kraus_gateindex=[]#List of Kraus index for each noise gate

        self._current_gateindex=0#Index of the current gate

    def add_x(self, qindex):
        self._exact_noiseless_circuit.x(qindex)
        self._exact_noisy_circuit.x(qindex)
        self._circuit_description.append(('x',qindex))
        self._current_gateindex+=1

    def add_y(self, qindex):
        self._exact_noiseless_circuit.y(qindex)
        self._exact_noisy_circuit.y(qindex)
        self._circuit_description.append(('y',qindex))
        self._current_gateindex+=1

    def add_z(self, qindex):
        self._exact_noiseless_circuit.z(qindex)
        self._exact_noisy_circuit.z(qindex)
        self._circuit_description.append(('z',qindex))
        self._current_gateindex+=1


    def add_s(self, qindex):
        self._exact_noiseless_circuit.s(qindex)
        self._exact_noisy_circuit.s(qindex)
        self._circuit_description.append(('s',qindex))
        self._current_gateindex+=1

    def add_cnot(self, qindex1, qindex2):
        self._exact_noiseless_circuit.cx(qindex1,qindex2)
        self._exact_noisy_circuit.cx(qindex1,qindex2)
        self._circuit_description.append(('cnot',(qindex1,qindex2)))
        self._current_gateindex+=1


    def add_hadamard(self, qindex):
        self._exact_noiseless_circuit.h(qindex)
        self._exact_noisy_circuit.h(qindex)
        self._circuit_description.append(('h',qindex))
        self._current_gateindex+=1


    #Inject a noise with Kgadget method
    def inject_noise(self, qindex):
        self._circuit_description.append(('K',qindex))
        self._Kraus_gateindex.append(self._current_gateindex)
        self._current_gateindex+=1



    #The state initilizarion circuit 
    def kadget_initialize_circuit(self,statestr=None):
        stateinit=QuantumCircuit(self._t,self._t)

        pass


    def compile_kgadget_circuit(self):
        dataqubit=QuantumRegister(self._n,name='data')
        noisequbit=QuantumRegister(self._t,name='noise')
        self._kgadget_circ= QuantumCircuit(dataqubit,noisequbit,)
        for i in range(0,self._t):

            self._kgadget_circ.append(Kraus(self._p),[noisequbit[i]])
            pass
        

    #Calculate the fidelity after run several simulations        
    def fidelity(self):
        pass


    def run(self):
        pass