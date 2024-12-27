# KGadget

Programmable K-gadget noisy quantum simulator





# Example of usage

I implement the KGadgetSimulator simulator class, where you can construct any Clifford circuit and specify some fixed location to add 
amplitude damping noise. The KgadgetSimulator class can calculate the exact noisy simulation, the exact KGadget simulation, and the compressed KGadget simulation. 



```python
ksim=KGadgetSimulator(1,1,0.5)
ksim.add_hadamard(0)
ksim.inject_noise(0)
ksim.add_z(0)
ksim.inject_noise(0)
ksim.compile_kgadget_circuit()

ksim.run_exact_kgadget_circuit()
ksim.run_compressed_kadget_circuit()
ksim.run_exact_noisy_circuit()

ksim.fidelity_of_exact_K_gadget()
```