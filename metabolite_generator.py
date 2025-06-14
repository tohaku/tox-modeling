try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError as exc:
    raise ImportError(
        "RDKit is required to run metabolite_generator. Please install the 'rdkit' package."
    ) from exc

# Simplified apply_reaction function (or copy the more robust one)
def apply_reaction(reaction_smarts, reactant_smiles_list):
    """
    Applies a reaction SMARTS to a list of reactant SMILES strings.
    Returns a list of product SMILES strings.
    """
    try:
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        
        # Create RDKit mol objects from SMILES and add explicit hydrogens
        # This can help with SMARTS patterns that explicitly match hydrogens (e.g., [c:1][H])
        reactant_mols = []
        for s in reactant_smiles_list:
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                # print(f"Warning: Invalid SMILES {s} for reaction {reaction_smarts}")
                return []
            reactant_mols.append(Chem.AddHs(mol))
        
        reactants_tuple = tuple(reactant_mols)
        products_sets = rxn.RunReactants(reactants_tuple)
        
        # Flatten the list of product tuples and convert to canonical SMILES
        # RemoveHs before generating SMILES for canonical representation if needed,
        # but MolToSmiles usually handles this well.
        # Using canonical SMILES helps in tracking unique molecules
        product_smiles_list = []
        for product_tuple in products_sets:
            for mol_in_tuple in product_tuple:
                if mol_in_tuple: # Ensure product molecule is valid
                    # Option: Remove Hs from product before canonical SMILES if AddHs was problematic for products
                    # mol_no_hs = Chem.RemoveHs(mol_in_tuple) 
                    # product_smiles_list.append(Chem.MolToSmiles(mol_no_hs, canonical=True))
                    product_smiles_list.append(Chem.MolToSmiles(mol_in_tuple, canonical=True))
        return product_smiles_list
    except Exception as e:
        # print(f"Error during reaction application: {e} with SMARTS {reaction_smarts} and reactants {reactant_smiles_list}")
        return []

REACTION_RULES = {
    'Aromatic_Hydroxylation': '[c:1][H:2]>>[c:1][O][H:2]',  # Matches an aromatic C-H bond
    'N_Dealkylation': '[N:1]([C:2][H])([*:3])([*:4])>>[N:1]([H])([*:3])([*:4]).[C:2]=O', # Simplified: N-CH becomes N-H + C=O
    'O_Glucuronidation': '[O:1][H]>>[O:1][C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@@H](C(=O)O)O1', # Attaches glucuronic acid to an OH group
    'Ester_Hydrolysis': '[C:1](=O)[O:2][C:3]>>[C:1](=O)[O:2][H].[C:3][O][H]', # Ester to acid and alcohol
}

def generate_metabolites_iterative(initial_smiles, reaction_rules_dict, iterations=1):
    """
    Generates metabolites by iteratively applying reaction rules.

    Args:
        initial_smiles (str): SMILES string of the initial molecule.
        reaction_rules_dict (dict): Dictionary of reaction_name: reaction_smarts.
        iterations (int): Number of iterations to apply reactions.

    Returns:
        set: A set of unique SMILES strings of all generated metabolites.
    """
    all_metabolites = {initial_smiles}
    current_generation_smiles = {initial_smiles}

    for i in range(iterations):
        next_generation_smiles = set()
        for smiles in current_generation_smiles:
            # Some reactions might require a co-reactant, but metabolic reactions are often enzymatic,
            # acting on the substrate directly or with implicit co-factors (like water, oxygen).
            # For rules needing H (e.g. Aromatic_Hydroxylation [c:1][H:2]), RDKit handles matching.
            # For rules producing two products (e.g. N_Dealkylation, Ester_Hydrolysis), both are collected.

            for rule_name, smarts in reaction_rules_dict.items():
                # Most metabolic reactions are unimolecular (substrate + enzyme/co-factor)
                # or pseudo-unimolecular (e.g., hydrolysis where water is abundant).
                # The provided SMARTS are written to take one reactant (the drug/metabolite).
                # If a rule implies a co-reactant (e.g. water for hydrolysis, oxygen for hydroxylation),
                # it's often implicit in the SMARTS transformation or handled by the enzyme's context.
                
                # Aromatic Hydroxylation: needs one reactant (the aromatic compound)
                # N-Dealkylation: needs one reactant (the N-alkylated compound)
                # O-Glucuronidation: needs one reactant (compound with -OH)
                # Ester Hydrolysis: needs one reactant (the ester). Water is implicit in the SMARTS products.

                products = apply_reaction(smarts, [smiles])
                for p_smiles in products:
                    if p_smiles not in all_metabolites:
                        next_generation_smiles.add(p_smiles)
                        all_metabolites.add(p_smiles)
        
        if not next_generation_smiles: # No new metabolites generated
            break
        current_generation_smiles = next_generation_smiles
        # print(f"Iteration {i+1}, found {len(next_generation_smiles)} new metabolites. Total unique: {len(all_metabolites)}")


    return all_metabolites


if __name__ == '__main__':
    # Test cases
    test_molecules = {
        "Toluene": "Cc1ccccc1",
        "Ethyl Acetate": "CC(=O)OCC",
        "Lidocaine-like": "CCN(CC)C(=O)c1c(C)cccc1C" # Simplified Lidocaine (no amide N-ethyl, different aromatic methyls)
    }

    iterations = 2 # Or 1, depending on desired depth

    for name, smiles in test_molecules.items():
        print(f"Generating metabolites for {name} ({smiles}) over {iterations} iterations:")
        
        # No need to pre-process with AddHs here if apply_reaction handles it.
        metabolites = generate_metabolites_iterative(smiles, REACTION_RULES, iterations)
        
        print(f"  Initial molecule: {smiles}")
        print(f"  Found {len(metabolites)} total unique molecules (initial + metabolites):")
        for i, met_smiles in enumerate(metabolites):
            print(f"    {i+1}. {met_smiles}")
        print("-" * 40)

    # Example of O-Glucuronidation on a molecule with an OH group (phenol from toluene hydroxylation)
    print("Specific test for O-Glucuronidation on a phenol:")
    phenol_smiles = "Oc1ccccc1" # A potential metabolite of benzene or toluene
    
    # apply_reaction now handles AddHs
    glucuronide_products = apply_reaction(REACTION_RULES['O_Glucuronidation'], [phenol_smiles])
    print(f"  Reactant (Phenol): {phenol_smiles}")
    print(f"  Glucuronide Products: {glucuronide_products}")
    print("-" * 40)
