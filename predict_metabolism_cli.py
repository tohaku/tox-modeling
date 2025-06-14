try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError as exc:
    raise ImportError(
        "RDKit is required to run predict_metabolism_cli. Please install the 'rdkit' package."
    ) from exc

REACTION_RULES = {
    'Aromatic_Hydroxylation': '[c:1][H:2]>>[c:1][O][H:2]',
    'N_Dealkylation': '[N:1]([C:2][H])([*:3])([*:4])>>[N:1]([H])([*:3])([*:4]).[C:2]=O',
    'O_Glucuronidation': '[O:1][H]>>[O:1][C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@@H](C(=O)O)O1',
    'Ester_Hydrolysis': '[C:1](=O)[O:2][C:3]>>[C:1](=O)[O:2][H].[C:3][O][H]',
}

def apply_reaction_for_iterative(reaction_smarts, reactant_smiles):
    """
    Applies a reaction SMARTS to a single reactant SMILES string.
    Adds explicit Hs to the reactant before reaction.
    Returns a list of product SMILES strings (canonical).
    """
    try:
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        mol = Chem.MolFromSmiles(reactant_smiles)
        if not mol:
            # print(f"Warning: Could not parse SMILES: {reactant_smiles}")
            return []

        mol_with_hs = Chem.AddHs(mol)
        products_sets = rxn.RunReactants((mol_with_hs,)) # Pass as a tuple of one reactant

        product_smiles_list = []
        for product_tuple in products_sets:
            for p_mol in product_tuple:
                if p_mol:
                    # Removing Hs from product before SMILES generation for cleaner output,
                    # but ensuring sanitization handles valencies correctly.
                    # Chem.SanitizeMol might be needed if RemoveHs causes issues.
                    # For now, let MolToSmiles handle it, as AddHs was on reactant.
                    # If product SMILES become problematic (e.g. retaining explicit Hs unnecessarily),
                    # consider Chem.RemoveHs(p_mol) here.
                    product_smiles_list.append(Chem.MolToSmiles(p_mol, canonical=True))
        return product_smiles_list
    except Exception as e:
        # print(f"Error in apply_reaction_for_iterative: {e}")
        return []

def generate_metabolites_iterative(initial_smiles, reaction_rules_dict, iterations=1):
    """
    Generates metabolites and records pathway steps.

    Args:
        initial_smiles (str): SMILES string of the initial molecule.
        reaction_rules_dict (dict): Dictionary of reaction_name: reaction_smarts.
        iterations (int): Number of iterations to apply reactions.

    Returns:
        tuple: (all_unique_smiles_set, pathway_steps_list)
               pathway_steps_list contains tuples of:
               (iteration, reactant_smiles, reaction_name, list_of_product_smiles)
    """
    all_unique_smiles = {initial_smiles}
    pathway_steps = []
    
    current_generation_smiles = {initial_smiles}

    for i in range(iterations):
        next_generation_smiles = set()
        if not current_generation_smiles: # Stop if no reactants for this generation
            break

        # print(f"Debug: Iteration {i+1}, processing {len(current_generation_smiles)} SMILES")
        for r_smiles in current_generation_smiles:
            for rule_name, smarts in reaction_rules_dict.items():
                # print(f"Debug: Applying {rule_name} to {r_smiles}")
                products = apply_reaction_for_iterative(smarts, r_smiles)
                # print(f"Debug: Products from {rule_name} on {r_smiles}: {products}")
                
                if products:
                    # Record the step, even if products are already seen,
                    # as this shows the transformation path.
                    pathway_steps.append((i + 1, r_smiles, rule_name, products))
                    
                    for p_smiles in products:
                        if p_smiles not in all_unique_smiles:
                            next_generation_smiles.add(p_smiles)
                            all_unique_smiles.add(p_smiles)
                        # If product is already in all_unique_smiles, but generated again
                        # from a different path or reactant, it could still be part of 
                        # next_generation_smiles if we want to explore its further reactions
                        # in this iteration block. However, current_generation_smiles is built
                        # from *newly added unique* smiles to avoid redundant processing loops
                        # on already processed structures within an iteration pass.
                        # The prompt implies we should process *all* products for next round.
                        # Let's ensure all products from this generation are considered for the next.
                        # The check `if p_smiles not in all_unique_smiles:` handles adding to the *overall* set.
                        # For `next_generation_smiles`, we add any product that could be a reactant.
                        # To avoid processing the same molecule multiple times *as a reactant* in the *next* iteration:
                        # We only add to next_generation_smiles if it's truly new, or if we want to allow
                        # a molecule to be processed again if it's formed via a new path.
                        # The current logic: next_generation_smiles gets only *newly discovered* unique SMILES.
                        # This means if X -> Y, and Z -> Y, Y is only added to next_generation_smiles once.
                        # This is generally what's desired to avoid exponential blowup on common intermediates.


        if not next_generation_smiles: # No new unique metabolites generated this round
            # print("Debug: No new unique metabolites generated this iteration.")
            break 
        current_generation_smiles = next_generation_smiles # Only process newly found unique smiles as reactants in next iter
        # print(f"Debug: End of Iteration {i+1}, next_generation_smiles for processing: {len(current_generation_smiles)}")

    return all_unique_smiles, pathway_steps

def main_cli():
    """
    Main command-line interface function.
    Prompts for SMILES and iterations, then prints metabolites and pathways.
    """
    print("Welcome to the Metabolism Predictor CLI!")
    
    initial_smiles = ""
    while not initial_smiles:
        initial_smiles_input = input("Enter the initial molecule SMILES string: ").strip()
        # Basic validation: check if RDKit can parse it
        mol = Chem.MolFromSmiles(initial_smiles_input)
        if mol:
            initial_smiles = initial_smiles_input # Use user's input style if valid
            # Or, canonicalize it early: initial_smiles = Chem.MolToSmiles(mol, canonical=True)
        else:
            print("Invalid SMILES string. Please try again.")

    num_iterations = -1
    while num_iterations < 0:
        try:
            num_iterations_input = input("Enter the number of metabolic generations (e.g., 1, 2, 3): ").strip()
            num_iterations = int(num_iterations_input)
            if num_iterations < 0:
                 print("Number of generations cannot be negative. Please enter a valid number.")
        except ValueError:
            print("Invalid input. Please enter an integer number.")

    print(f"\nProcessing {initial_smiles} for {num_iterations} generation(s)...")

    all_metabolites_set, pathway_details_list = generate_metabolites_iterative(
        initial_smiles, 
        REACTION_RULES, 
        num_iterations
    )

    print("\n--- All Unique Molecules (Initial + Metabolites) ---")
    if all_metabolites_set:
        for idx, smiles in enumerate(sorted(list(all_metabolites_set))): # Sorted for consistent output
            print(f"{idx + 1}. {smiles}")
    else:
        print("No metabolites generated or initial SMILES was invalid.")

    print("\n--- Metabolic Pathway Steps ---")
    if pathway_details_list:
        for step_num, (iteration, r_smiles, rule_name, p_smiles_list) in enumerate(pathway_details_list):
            print(f"\nStep {step_num + 1} (Generation {iteration}):")
            print(f"  Reactant: {r_smiles}")
            print(f"  Reaction: {rule_name}")
            print(f"  Products:")
            for p_idx, p_smiles in enumerate(p_smiles_list):
                print(f"    {p_idx + 1}. {p_smiles}")
    else:
        print("No metabolic steps were generated (or 0 iterations selected).")
    
    print("\nMetabolism prediction finished.")

if __name__ == '__main__':
    main_cli()
