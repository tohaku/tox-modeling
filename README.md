# Tox Modeling

This repository contains small RDKit-based utilities for exploring basic chemical reaction predictions.

- **metabolite_generator.py** iteratively applies simple metabolic reaction rules to a starting molecule and prints the resulting metabolites.
- **reaction_mapper.py** demonstrates applying named reaction SMARTS to reactants for educational purposes.
- **predict_metabolism_cli.py** provides a small command line interface for generating metabolites and tracing pathways.

RDKit must be installed for these scripts to run. The provided `requirements.txt` lists the dependency.

## Usage

1. Install dependencies:

   ```bash
   pip install -r requirements.txt
   ```

2. Run the example metabolite generator:

   ```bash
   python metabolite_generator.py
   ```

   This will print sample metabolites for a few hard-coded molecules.

3. Map reactions using the example rules:

   ```bash
   python reaction_mapper.py
   ```

   The script applies a handful of named reactions to sample reactants and prints the products.

4. Start the interactive metabolism prediction CLI:

   ```bash
   python predict_metabolism_cli.py
   ```

   You will be prompted for an initial SMILES string and the number of generations to simulate. The CLI displays the unique products and the reaction pathway.

