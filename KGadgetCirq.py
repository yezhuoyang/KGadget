import cirq
import numpy as np
from qiskit.quantum_info import random_statevector, Statevector
from qiskit.quantum_info import random_clifford
from qiskit.quantum_info.operators import Operator
import json

class Proj(cirq.Gate):
    def __init__(self, on_0: bool=True):
        # on_0 = True means project to 0. If False, project to 1.
        self.on_0 = on_0
        super(Proj, self)

    def _num_qubits_(self):
        return 1

    def _unitary_(self):
        # Doesn't have to be unitary, this is just the matrix
        # that will be applied at runtime!
        return np.diag([1,0]) if self.on_0 else np.diag([0,1])

    def _circuit_diagram_info_(self, args):
        return "P0"


class RandomCliff(cirq.Gate):
    def __init__(self,n):
        # on_0 = True means project to 0. If False, project to 1.
        self.n = n
        self.data = Operator(random_clifford(num_qubits=self.n)).data
        super(RandomCliff, self)

    def _num_qubits_(self):
        return self.n

    def _unitary_(self):
        # Doesn't have to be unitary, this is just the matrix
        # that will be applied at runtime!
        return self.data

    def _circuit_diagram_info_(self, args):
        return ("C",) * self.n


class KGadgetSimulator:
    def __init__(self,n,t,p):
        #Define Paramters
        self.n = n
        self.t = t
        self.p = p

        #Instantiate Projection and K
        self.proj = Proj(on_0=True)
        self.initial_psi = [random_statevector(2).data for _ in range(n)] #We construct the state vector by storing individual qubit data, since storing 2**n vector is too expensive
        self.initial_K = np.array([1, np.sqrt(1 - self.p)])
        self.initial_K /= np.linalg.norm(self.initial_K)  # Ensure normalization
        
        #Generate a Random Clifford Sequence
        self.clifford_seq = [RandomCliff(n) for _ in range(t)]


    #GENERAL HELPERS
    @staticmethod
    def kron(args):
        current = args[0]
        for s in args[1:]:
            current = np.kron(current, s)
        return current


    #COMPRESSION HELPERS
    @staticmethod
    def standard_basis(n):
        return [list(np.eye(n, dtype=int)[i]) for i in range(n)] + [[0 for _ in range(n)]]

    @staticmethod
    def bit_conversion(i):
        if i==1:
            return np.array([1,0])
        elif i==0:
            return np.sqrt(1/2)*np.array([1,1])


    #SIMULATORS
    def exact_circuit(self,):
        '''
        returns circuit of the exact K-Gadget Simulation
        '''

        # Define the qubits
        psi_qubits = [cirq.LineQubit(i) for i in range(n)]
        K_qubits = [cirq.LineQubit(i+n) for i in range(t)] 

        # Create a quantum circuit
        circuit = cirq.Circuit()

        # Instantiate K state for each K qubit
        for K in K_qubits:
            circuit.append(cirq.StatePreparationChannel(self.initial_K).on(K))

        # Instantiate a random psi state
        for psi,psi_data in zip(psi_qubits,self.initial_psi):
            circuit.append(cirq.StatePreparationChannel(psi_data).on(psi))

        #Apply Gadgets
        for K,cliff in zip(K_qubits,self.clifford_seq):
            circuit.append(cliff.on(*psi_qubits)) #apply intermediary clifford on psi
            circuit.append(cirq.CNOT(psi_qubits[0], K))  #assuming that K is applied to the first qubit (psi0)
            circuit.append(self.proj.on(K))
            circuit.append(cirq.measure(K, key='m'))

        return circuit


    def compressed_partial_circuit(self, bit_str):
        '''
        returns circuit of the basis compressed K-Gadget Simulation
        '''

        # Define the qubits
        psi_qubits = [cirq.LineQubit(i) for i in range(n)]
        K_qubits = [cirq.LineQubit(i+n) for i in range(t)] 

        # Create a quantum circuit
        circuit = cirq.Circuit()

        #Prepare the partial compressed state
        for K,bit in zip(K_qubits,bit_str):
            vector = self.bit_conversion(bit)
            circuit.append(cirq.StatePreparationChannel(vector).on(K))

        # Instantiate a random psi state
        for psi,psi_data in zip(psi_qubits,self.initial_psi):
            circuit.append(cirq.StatePreparationChannel(psi_data).on(psi))


        #Apply Gadgets
        for K,cliff in zip(K_qubits,self.clifford_seq):
            circuit.append(cliff.on(*psi_qubits)) #apply intermediary clifford on psi
            circuit.append(cirq.CNOT(psi_qubits[0], K))  #assuming that K is applied to the first qubit (psi0)
            circuit.append(self.proj.on(K))
            circuit.append(cirq.measure(K, key='m'))

        return circuit


    #SIMULATORS 
    def compressed_simulation(self,):
        alpha = (1-np.sqrt(1-self.p))/(np.sqrt(2*(1-self.p)))
        basis = self.standard_basis(self.t)

        stab_states = []
        for bit_str in basis:
            alpha_scalar = alpha**(np.sum(bit_str))
            circuit = self.compressed_partial_circuit(bit_str)
            raw_stab_state = self.run_simulation(circuit)
            stab_states.append(raw_stab_state*alpha_scalar)

        output_state = np.sum(np.array(stab_states),axis=0)
        output_state = output_state/np.linalg.norm(output_state) #norm final output state
        return output_state


    def run_simulation(self,circuit):
        simulator = cirq.Simulator(dtype=np.complex128)
        result = simulator.simulate(circuit)

        #Extract Simulated Ouput State
        indices_to_keep = list(range(0, 2**(self.n+self.t), 2**self.t))
        output_state = result.final_state_vector[indices_to_keep]
        return output_state

    def true_sim(self,):
        '''
        return exact state vector simulation
        '''

        true = self.kron(self.initial_psi)
        for cliff in self.clifford_seq:
            K_n = [np.diag(self.initial_K)] + [np.eye(2) for _ in range(n-1)] #K is only applied on the first qubit
            true = cliff.data@true
            true = self.kron(K_n)@true
        true = true/np.linalg.norm(true)
        return true


    #ERROR STATISTICS
    def get_error(self,x,y):
        return 1 - abs(np.vdot(x,y))**2



def plot_errors(errors,y,log=False):

    data = list(errors.values())

    averages = [np.mean(lst) for lst in data]
    mins = [np.min(lst) for lst in data]
    maxs = [np.max(lst) for lst in data]

    lower_error = [avg - mn for avg, mn in zip(averages, mins)]
    upper_error = [mx - avg for avg, mx in zip(averages, maxs)]
    errors = [lower_error, upper_error]

    if log:
        errors = [np.abs(np.log(lower_error)), np.abs(np.log(upper_error))]
        print(min(np.log(lower_error)))
        print(min(np.log(upper_error)))
        plt.errorbar(range(1,len(data)+1), np.log(averages), yerr=errors, fmt='o', capsize=5, capthick=2, label='Average with Min-Max Bars')
        plt.plot(range(1,len(data)+1), np.log(averages), linestyle='-', color='b')  # 'b' specifies the color blue; adjust as needed
    else:
        plt.errorbar(range(1,len(data)+1), averages, yerr=errors, fmt='o', capsize=5, capthick=2, label='Average with Min-Max Bars')
        plt.plot(range(1,len(data)+1), averages, linestyle='-', color='b')  # 'b' specifies the color blue; adjust as needed


    plt.xlabel(y)
    plt.ylabel('error')
    plt.show()


def K_Fidelty(t,a,v):
    #Calculate error for approximations made using the basis vectors and 0
    #The Proof of this can be found in Circuit Depth / Size Investigation
    K_L =  2*t*a*v + t*(a**2) + (t**2-t)*(a**2)*(v**2) + 1
    K_F = (a**2 + 2*a*v + 1)**t
    N = (K_L*K_F)**(-1/2)
    S = (1+a*v)**t + t*(a**2 + a*v)*(1+a*v)**(t-1)
    delta = 1 - (S*N)**2
    return delta




if __name__ == '__main__':
    import matplotlib.pyplot as plt 
    from tqdm import tqdm
    
    '''
    #SIMULATIONS WITH RESPECT TO n AND t
    N = 100
    max_t = 21
    n = 1
    p = .10

    errors_t = {i:[] for i in range(1,max_t)}
    for _ in tqdm(range(N)):
        for t in range(1,max_t):
            sim = KGadgetSimulator(n,t,p)
            compressed_output = sim.compressed_simulation()
            true_output = sim.true_sim()
            delta = sim.get_error(compressed_output,true_output)
            errors_t[t].append(float(delta))

    with open('errors_t.json', 'w') as json_file:
        json.dump(errors_t, json_file, indent=4)
    
        
    N = 100
    max_n = 11
    t = 10
    p = .10
    errors_n = {i:[] for i in range(1,max_n)}
    for _ in tqdm(range(N)):
        for n in tqdm(range(1,max_n)):
            sim = KGadgetSimulator(n,t,p)
            compressed_output = sim.compressed_simulation()
            exact_circuit = sim.exact_circuit()
            exact_output = sim.run_simulation(exact_circuit)

            delta = sim.get_error(compressed_output,exact_output)
            errors_n[n].append(float(delta))

    with open('errors_n.json', 'w') as json_file:
        json.dump(errors_n, json_file, indent=4)
    
    
    

    # PLOTTING 
    with open('errors_t.json', 'r') as json_file:
        errors = json.load(json_file)

    plot_errors(errors,'t',log=False)
    plot_errors(errors,'t',log=True)

    
    with open('errors_n.json', 'r') as json_file:
        errors = json.load(json_file)
        
    plot_errors(errors,'n')
    plot_errors(errors,'t',log=True)
    



    #PLOTTING ERROR WITH RESPECT TO FIDELITY 
    with open('errors_t.json', 'r') as json_file:
        errors = json.load(json_file)

    fids = []
    errors_avgs = []
    for t in range(1,len(errors)):
        p = .10
        alpha = (1-np.sqrt(1-p))/(np.sqrt(2*(1-p)))
        v = np.sqrt(1/2)
        fidelity = K_Fidelty(t,alpha,v)

        fids.append(fidelity)
        errors_avgs.append(np.mean(errors[str(t)]))

    plt.plot(fids,errors_avgs)
    plt.xlabel('K Fidelity')
    plt.ylabel('Error')
    plt.show()
    '''
    

    
    #EXAMPLE USAGE
    p = 0.10
    n = 1  # Number of qubits for psi
    t = 3  # Number of K qubits 

    sim = KGadgetSimulator(n,t,p)

    #Compressed Simulation
    compressed_output = sim.compressed_simulation()

    #Exact Simulation
    exact_circuit = sim.exact_circuit()
    exact_output = sim.run_simulation(exact_circuit)

    print(exact_circuit)

    #True Simulation
    #true_output = sim.true_sim()

    #Error statistics
    exact_error = sim.get_error(exact_output,compressed_output)
    print(exact_error)
    