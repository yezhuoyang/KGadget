from qiskit import QuantumCircuit, transpile, QuantumRegister, ClassicalRegister
from qiskit_aer import AerSimulator
from  qiskit.quantum_info import Kraus
from qiskit.quantum_info import Statevector
from qiskit.quantum_info import Operator
from qiskit.quantum_info import partial_trace
import math


#Label convention for K gadget simulator:
# Qubit 0~n-1: Data qubits, Qubit n~n+t-1: Noise qubits
# But keep in mind that the qiskit index is reversed, so the qubit 0 is the rightmost qubit in the circuit
# The state vector is |n_t-1,n_t-2,...,n_0,q_n-1,q_n-2,...,q_0>
# To act a singe qubit gate on the kth qubit, the operator is   I^{n-k-1} ^ op ^ I^{k}, note that the order of the operator is is reversed
# If we use statevector operator string, the string shoould be 'I'*(n-k-1) + op + 'I'*k
# To act two qubits CNOT gate on qubits (k1,k2), the operator is a two qubit gate with CNOT^{k1,k2}


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



proj1All = lambda n,t: Operator.from_label('1'*t)^iexpand(n)



#Act gate on a qubit
single_op_on_qubit = lambda opstr,ntotal,qubitindex: Operator.from_label('I'*(ntotal-qubitindex-1) + opstr +  'I'*(qubitindex))


proj_on_qubit = lambda projop,ntotal,qubitindex: (
                projop if ntotal==1 else
                iexpand(ntotal-1)^projop  if qubitindex == 0 else
                projop^iexpand(qubitindex) if qubitindex == ntotal - 1 else
                iexpand(ntotal-qubitindex-1)^ projop ^ iexpand(qubitindex)
)


Kop_onqubit = lambda p, ntotal, qubitindex: (
    Kop(p) if ntotal==1 else
    iexpand(ntotal - qubitindex - 1)^Kop(p)   if qubitindex == 0 else
    Kop(p)^iexpand(qubitindex)   if qubitindex == ntotal - 1 else
     iexpand(ntotal - qubitindex - 1)^ Kop(p) ^iexpand(qubitindex) 
)


def cnot_on_qubits(ntotal,qindex1,qindex2):
    circ=QuantumCircuit(ntotal)
    circ.cx(qindex1,qindex2)
    return Operator.from_circuit(circ)


def normalize(statevector):
    norm=abs(statevector.inner(statevector))**0.5
    return statevector/norm


def fidelity(statevector1,statevector2):
    return abs(statevector1.inner(statevector2))**2



def theoretical_fidelity(alpha, s, p):
    """
    Computes the function:
    
        Fid(alpha)^(s,p) = (1 / (1 + p)^s)
                           * (2^(-s))
                           * (1 / (1 + [s*(s+1)/2]*alpha^2 + sqrt(2)*s*alpha))
                           * { ( sqrt(2)*alpha*s + sqrt(p) + 1 ) * ( sqrt(p)+1 )^(s-1) }^2
    
    Parameters
    ----------
    alpha : float
        The alpha parameter (real).
    s : float
        The s parameter (real).
    p : float
        The p parameter (real).
    
    Returns
    -------
    float
        The value of the function for the given alpha, s, and p.
    """
    # Term 1: (1 / (1+p)^s)
    term1 = 1 / ((1 + p) ** s)
    
    # Term 2: 2^(-s)
    term2 = 2 ** (-s)
    
    # Denominator in fraction: (1 + s(s+1)/2 * alpha^2 + sqrt(2)*s*alpha)
    denom = 1 + (s * (s + 1) / 2) * (alpha ** 2) + math.sqrt(2) * s * alpha
    
    # Term 3: 1 / denom
    term3 = 1 / denom
    
    # Inner bracket: ( sqrt(2)*alpha*s + sqrt(p) + 1 ) * ( sqrt(p)+1 )^(s-1)
    bracket = (math.sqrt(2) * alpha * s + math.sqrt(p) + 1) * ((math.sqrt(p) + 1) ** (s - 1))
    
    # Term 4: [ bracket ]^2
    term4 = bracket ** 2
    
    return term1 * term2 * term3 * term4





#Qiskit doesn't not support non-unitary circuit simulation
#Noisy Quantum Circuit which simulate a quantum circuit with non-unitary noise
#by naive state vector simulation
class noisyQcircuit:

    def __init__(self,n,t,p):
        self._nqubit=n
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


def reverse_binary_string(x: int) -> str:
    """
    Convert a binary integer x to the reversed version of its binary string 
    (excluding the '0b' prefix), preserving the original bit length.
    
    Examples:
    ----------
    reverse_binary_string(0b1110000)   -> '0000111'
    reverse_binary_string(0b10100100)  -> '00100101'
    """
    # Convert integer x to a binary string (e.g., '0b1110000' -> '1110000')
    bin_str = bin(x)[2:]
    # Return the reversed binary string
    return bin_str[::-1]


class KGadgetSimulator:


    def __init__(self,n,t,p):
        #Define Paramters
        self._n = n#Number of qubits
        self._t = t#Number of different noise
        self._p = p#Parameter of amplititude damping noise

        self._exact_noiseless_simulator = AerSimulator(method="stabilizer")#Qiskit stabilizer simulator for exact simulation without any noise
        self._exact_noiseless_circuit=QuantumCircuit(n,n)


        self._exact_noise_simulator = noisyQcircuit(n,t,p)#Qiskit stabilizer simulator for exact simulation with  noise
        self._exact_noisy_final_state=None#Final state vector of exact simulation with noise


        self._kgadget_circ=None#Circuit for Kgadget simulation
        self._kgadget_simulator= AerSimulator(method="stabilizer")#Qiskit stabilizer simulator for Kgadget simulation
        self._kgadget_compressed_final_state=None
        self._kgadget_non_compressed_final_state=None


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
        kstate=Statevector([1,self._p**0.5])
        result=kstate.copy()
        for i in range(self._t-1):
            result=result.tensor(kstate)

        zeros=Statevector.from_label('0'*self._n)
        result=result.tensor(zeros)    
        return normalize(result)    #Normalize the state vector



    def compile_kgadget_circuit(self,initState=None):
        dataqubit=QuantumRegister(self._n,name='data')
        noisequbit=QuantumRegister(self._t,name='noise')
        self._kgadget_circ= QuantumCircuit(dataqubit,noisequbit)
        if initState is not None:
            self._kgadget_circ.set_statevector(initState)
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


    #Run the exact Kgdaget scheme
    def run_exact_kgadget_circuit(self):
        initial_state=self.exact_kstates()

        self.compile_kgadget_circuit(initial_state)
        self._kgadget_circ.save_statevector()

        stabilizer_simulator = AerSimulator(method="statevector")

        circ = transpile(self._kgadget_circ, stabilizer_simulator)
        # Run and get statevector
        result = stabilizer_simulator.run(circ).result()
        statevector = result.get_statevector(self._kgadget_circ)

        # Project all helper qubits to |1>
        statevector=statevector.evolve(proj1All(self._n,self._t))


        self._kgadget_non_compressed_final_state=normalize(statevector)
        return self._kgadget_non_compressed_final_state





    #Run the compressed Kgdaget scheme with a given input string
    def run_compressed_kadget_circuit_single(self,inputstr):
        simcirc=QuantumCircuit(self._n+self._t)
        for i in range(self._t):
            if inputstr[i]=='+':
                simcirc.h(self._n+i)
        self.compile_kgadget_circuit()
        simcirc.append(self._kgadget_circ,range(self._n+self._t))
        simcirc.save_statevector()

        stabilizer_simulator = AerSimulator(method="statevector")
        circ = transpile(simcirc, stabilizer_simulator)
        # Run and get statevector
        result = stabilizer_simulator.run(circ).result()
        statevector = result.get_statevector(simcirc)

        # Project all helper qubits to |1>
        statevector=statevector.evolve(proj1All(self._n,self._t))

        #print(statevector.to_dict())
        return statevector
        


    def run_compressed_kadget_circuit(self,alpha=1):
        finalstate=self.run_compressed_kadget_circuit_single('+'*self._t)
        for i in range(self._t):
            inputstr='0'*i+'+'+'0'*(self._t-i-1)
            finalstate+=alpha*self.run_compressed_kadget_circuit_single(inputstr)
        #normalize final state
        self._kgadget_compressed_final_state=normalize(finalstate)
        return self._kgadget_compressed_final_state

        

    #Calculate the fidelity after run several simulations  and take the average of the state vector     
    def fidelity_of_compression(self):
        return fidelity(self._kgadget_non_compressed_final_state,self._kgadget_compressed_final_state)




    #Calculate the fidelity of the exact K gadget method
    def fidelity_of_exact_K_gadget(self):
        kgadstate=self._kgadget_non_compressed_final_state.to_dict()
        print(kgadstate)
        exactstate=self._exact_noisy_final_state.to_dict()
        print(exactstate)
        fidelity=0
        for i in range(2**self._n):
            fidelity+=kgadstate['1'*self._t+reverse_binary_string(i)]*exactstate[reverse_binary_string(i)]
        return abs(fidelity)**2




    #Run the exact circuit without any noise and return a final state vector
    def run_exact_noiseless(self):
        self._exact_noiseless_circuit.save_statevector()

        stabilizer_simulator = AerSimulator(method="statevector")
        circ = transpile(self._exact_noiseless_circuit, stabilizer_simulator)
        # Run and get statevector
        result = stabilizer_simulator.run(circ).result()
        statevector = result.get_statevector(self._exact_noiseless_circui)

        self._exact_noiseless_final_state=normalize(statevector)
        return self._exact_noiseless_final_state


        


    def run_exact_noisy_circuit(self):
        self._exact_noisy_final_state=normalize(self._exact_noise_simulator.run())
        return self._exact_noisy_final_state
    





class repetitionQEC:


    def __init__(self,d,t,p):    
        self.kgadgetsim=KGadgetSimulator(d,t,p)




class surfaceQEC:


    def __init__(self,d,t,p):
        self.kgadgetsim=KGadgetSimulator(d,t,p)





if __name__ == '__main__':
    ksim=KGadgetSimulator(1,2,0.5)
    ksim.add_hadamard(0)
    ksim.inject_noise(0)
    ksim.add_z(0)
    ksim.inject_noise(0)

    #ksim.compile_kgadget_circuit()
    ksim.run_exact_kgadget_circuit()
    ksim.run_compressed_kadget_circuit()

    ksim.run_exact_noisy_circuit()

    ksim.fidelity_of_exact_K_gadget()


    #print(ksim.fidelity_of_compression())
    #print(ksim.run_kadget_circuit('++'))