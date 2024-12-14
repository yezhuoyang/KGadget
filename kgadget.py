from qiskit import QuantumCircuit, transpile, QuantumRegister, ClassicalRegister
from qiskit_aer import AerSimulator



class KGadgetSimulator:


    def __init__(self,n,t,p):
        #Define Paramters
        self.n = n#Number of qubits
        self.t = t#Number of different noise
        self.p = p#Parameter of amplititude damping noise
        self.simulator = AerSimulator(method="stabilizer")

        self.dataqubit=QuantumRegister(n,name='data')
        self.noisequbit=QuantumRegister(t,name='noise')

        self.circ= QuantumCircuit(self.dataqubit,self.noisequbit)

    #Inject a noise with Kgadget method
    def inject_noise(self, qindex):
        pass 


    def compile_circuit(self):
        pass



    def run(self):

        circ = transpile(self.circ, self.simulator)

        # Run and get statevector
        result = self.simulator.run(circ).result()

        return result