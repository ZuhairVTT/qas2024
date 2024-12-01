# Binding Energy Analysis

This notebook analyzes the binding energy from hybrid quantum-classical calculations by processing CP2K output files.

## Implementation

```python
import os
import re

def extract_energy_from_cp2k_log(filepath):
    """Extract the total energy from a CP2K log file using the FORCE_EVAL energy pattern."""
    try:
        with open(filepath, 'r') as f:
            content = f.read()
            # Look for the FORCE_EVAL total energy
            pattern = r'ENERGY\|\s*Total FORCE_EVAL \( QS \) energy \[a\.u\.\]:\s*(-?\d+\.\d+)'
            match = re.search(pattern, content)
            if match:
                return float(match.group(1))
            raise ValueError(f"Could not find FORCE_EVAL total energy in {filepath}")
    except FileNotFoundError:
        print(f"File not found: {filepath}")
        return None    

# Define paths to the log files
supercell_log = "2_supercell/cp2k.log"
substrate_log = "3_Al/cp2k.log"
inhibitor_log = "4_inhibitor/cp2k.log"

# Extract energies
E_supercell = extract_energy_from_cp2k_log(supercell_log)
E_substrate = extract_energy_from_cp2k_log(substrate_log)
E_inhibitor = extract_energy_from_cp2k_log(inhibitor_log)

# Calculate binding energy
if all(e is not None for e in [E_supercell, E_substrate, E_inhibitor]):
    E_binding = E_supercell - (E_substrate + E_inhibitor)
    
    print(f"Energies (Hartree):")
    print(f"Supercell: {E_supercell:.6f}")
    print(f"Substrate: {E_substrate:.6f}")
    print(f"Inhibitor: {E_inhibitor:.6f}")
    print(f"\nBinding Energy:")
    print(f"E_binding = {E_binding:.6f} Hartree")
    print(f"E_binding = {E_binding * 27.211386245988:.6f} eV")
    print(f"E_binding = {E_binding * 627.509474:.6f} kcal/mol")
```

## Results

The analysis produces the following results:

```
Energies (Hartree):
Supercell: -241.400215
Substrate: -198.415358
Inhibitor: -42.970690

Binding Energy:
E_binding = -0.014167 Hartree
E_binding = -0.385508 eV
E_binding = -8.890018 kcal/mol
```

## Key Components

1. **Energy Extraction**: Uses regular expressions to parse CP2K log files and extract total energies
2. **Unit Conversion**: Converts binding energy to different units:
   - Hartree (atomic units)
   - Electron volts (eV)
   - Kilocalories per mole (kcal/mol)
3. **Error Handling**: Includes checks for file existence and proper energy extraction

## Analysis Details

The binding energy calculation follows these steps:

1. Extract total energies from CP2K output files for:
   - Complete system (supercell)
   - Substrate alone
   - Inhibitor molecule alone
   
2. Calculate binding energy using the formula:
   ```
   E_binding = E_supercell - (E_substrate + E_inhibitor)
   ```

3. Convert the binding energy to different units using conversion factors:
   - 1 Hartree = 27.211386245988 eV
   - 1 Hartree = 627.509474 kcal/mol 