---
title: Quantum Computing for Materials Science
summary: Introduction to simulating materials using quantum computers
---

# Quantum Computing for Materials Science

## Contents

- [Introduction](#introduction)
- [Why Quantum Computing for Materials?](#why-quantum-computing-for-materials)
- [Hands-on Examples](#hands-on-examples)

## Introduction

Quantum chemistry calculations are among the most promising early applications of quantum computers. However, materials science simulations, particularly those involving solid-state systems, have not been well covered in the literature. Such calculations are usually associated with periodic system descriptions and can involve a large number of atoms, often in the hundreds. This repository aims to provide a comprehensive guide on how to perform these calculations using quantum computers in a way that aligns with current efforts to introduce quantum-centric supercomputers, where the quantum computer is considered an accelerator for classical calculations, similar to the role of GPUs in various fields. The tutorial here is focusing on simulating the calculation on a quantum computer, some hinders weren't solved when it comes to running the calculation on quantum computer at the time of writing this tutorial.

## Why Quantum Computing for Materials?

The simulation of large solid-state systems typically involves a significant number of atoms, which necessitates approximations in calculating the ground state energy. Methods like Density Functional Theory (DFT) are useful in this context. However, there are important parts of the system that can be described more accurately, especially when they involve interactions of particular interest to scientists. Quantum computing algorithms could be advantageous in focusing on such subsystems of the supercell and communicating the results to the DFT code. In simpler terms, they can correct the energy by accounting for more detailed information from these specific subsystems. This approach is known as quantum embedding[^1]. Of course, this is a very simplified description and can be considered analogous to the QM/MM method; again, this is just a simplification of the problem description and not an entirely accurate illustration. 

It is a hot topic now given the quantum centeric supercomputers effort and the rise of hybrid quantum classical approaches. More recent literature to study can be found in those references[^2] [^3] [^4] [^5].

## Hands-on Example

### 1. Classical calculations

The classical calculations are organized in several steps:

#### a. Supercell generation
- Setting up the supercell that contains the periodic system ready for the CP2K calculations. Example can be a supercell of Al(111) surface and inhibitor on top of the Al(111) substrate
- The script for generating the supercell is [`here`](0_classical_calculations/0_supercell_generation/supercell_generator) (download [`here`](0_classical_calculations/0_supercell_generation/supercell_generator.py)).

#### b. Geometry Optimization
- We need to optimize the geometry of the system to get the best possible structure to feed to the DFT calculation. Here we use ASE to optimize the structure with machine learning potential to accelerate the calculation with near DFT quality.
- More details regarding the geometry optimization steps are [`detailed here`](0_classical_calculations/1_geometry_optimization/gemetry_optimize.md) and the related script can be downloaded from [`here`](0_classical_calculations/1_geometry_optimization/geo_opt.py)

#### c. Supercell Calculations
- Now, we can perform the DFT calculations for the periodic system. Please follow the details [`here`](0_classical_calculations/2_supercell_calculation/supercell_calculation.md) regarding the system setup and calculational details where the supercell DFT calculations is performed using CP2K code.

#### d. Binding Energy Calculation
- Now, we can calculate the binding energy of the inhibitor on the Al(111) surface. The script for calculating the binding energy and more details are [`in this page`](0_classical_calculations/5_binding_energy/binding_energy.md)


### 3. Quantum System Preparation
- Example: [`2_supercell/`](2_supercell)
  ```python
  # Supercell creation for periodic boundary conditions
  from ase.build import make_supercell
  supercell = make_supercell(atoms, [[2,0,0], [0,2,0], [0,0,1]])
  ```

### 4. Quantum Simulations
- Al surface modeling ([`3_Al/`](3_Al))
- Inhibitor interaction studies ([`4_inhibitor/`](4_inhibitor))
- Example from [`analysis.ipynb`](analysis):
  ```python
  # Example adaptive VQE setup
  from qiskit.algorithms.minimum_eigensolvers import VQE
  from qiskit.circuit.library import TwoLocal
  
  ansatz = TwoLocal(num_qubits, 'ry', 'cz')
  vqe = VQE(ansatz, optimizer)
  ```

### 5. Practical Implementation Steps
1. System preparation
   - Geometry optimization
   - Electronic structure calculation
2. Quantum circuit construction
   - Hamiltonian mapping
   - Ansatz selection
3. VQE execution
   - Parameter optimization
   - Energy calculation
4. Results analysis
   - Convergence checking
   - Property calculation

### 7. Performance Comparisons
- Classical DFT vs Quantum results
- Resource requirements analysis
- Error analysis and mitigation strategies


## References

[^1]: Vorwerk, C., Sheng, N., Govoni, M., et al. "Quantum embedding theories to simulate condensed systems on quantum computers." *Nature Computational Science* **2**, 424–432 (2022). [DOI:10.1038/s43588-022-00279-0](https://doi.org/10.1038/s43588-022-00279-0)

[^2]: Multiscale Embedding for Quantum Computing. *arXiv preprint* arXiv:2409.06813 (2023). [https://arxiv.org/abs/2409.06813](https://arxiv.org/abs/2409.06813)

[^3]: Toward practical quantum embedding simulation of realistic chemical systems on near-term quantum computers. *Chemical Science* **13**, 7311-7324 (2022). [DOI:10.1039/D2SC01492K](https://pubs.rsc.org/en/content/articlelanding/2022/sc/d2sc01492k)

[^4]: Quantum Embedding and Quantum Simulations. Galli Group Research. University of Chicago. [https://galligroup.uchicago.edu/Research/embedding.php](https://galligroup.uchicago.edu/Research/embedding.php)

[^5]: Quantum Embedding Method for the Simulation of Strongly Correlated Systems on Quantum Computers. *arXiv preprint* arXiv:2302.03052 (2023). [https://arxiv.org/abs/2302.03052](https://arxiv.org/abs/2302.03052)