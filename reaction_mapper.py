from rdkit import Chem
from rdkit.Chem import AllChem

REACTION_RULES = {
    'Amide_Formation': '[C:1](=O)[OH1:2].[N;H2;!$(N=[*]):3]>>[C:1](=O)[N;H1:3].[OH2:2]', # Carboxylic acid + Primary Amine -> Amide + Water
    'Esterification': '[C:1](=O)[OH1:2].[C:3][OH1:4]>>[C:1](=O)[O:4][C:3].[OH2:2]', # Carboxylic acid + Alcohol -> Ester + Water
    'Friedel_Crafts_Acylation': '[c:1]Cl.[c:2](=O)Cl>>[c:1][c:2](=O).[Cl]', # This rule is problematic, see notes
    'Grignard_Reaction': '[C:1](=O).[C:2][MgX]>>[C:1](O[MgX])[C:2]', # Aldehyde/Ketone + Grignard
    'Diels_Alder': '[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2]=[C:3][C:4][C:5][C:6]1',
    'Williamson_Ether_Synthesis': '[C:1][O-].[C:2][Br]>>[C:1]O[C:2].[Br-]'
}

def apply_reaction(reaction_name, reactants_smiles):
    """
    Applies a named reaction to a list of reactant SMILES strings.

    Args:
        reaction_name (str): The name of the reaction from REACTION_RULES.
        reactants_smiles (list): A list of SMILES strings for the reactants.

    Returns:
        list: A list of SMILES strings for the products, or None if reaction fails.
    """
    if reaction_name not in REACTION_RULES:
        print(f"Error: Reaction '{reaction_name}' not found in REACTION_RULES.")
        return None

    reaction_smarts = REACTION_RULES[reaction_name]
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)

    reactants = [Chem.MolFromSmiles(s) for s in reactants_smiles]
    if None in reactants:
        print("Error: Invalid SMILES string provided for one or more reactants.")
        return None

    products_sets = rxn.RunReactants(tuple(reactants))

    if not products_sets:
        print(f"Reaction '{reaction_name}' did not yield any products with the given reactants.")
        return []

    # Flatten the list of product tuples and convert to SMILES
    product_smiles_list = []
    for product_tuple in products_sets:
        for mol in product_tuple:
            product_smiles_list.append(Chem.MolToSmiles(mol))
            
    return product_smiles_list

if __name__ == '__main__':
    test_cases = {
        "Benzene": ("Friedel_Crafts_Acylation", ["c1ccccc1Cl", "CC(=O)Cl"]), # Incorrect reactant for typical Friedel-Crafts
        "Trimethylamine": ("Amide_Formation", ["CC(=O)O", "CN(C)C"]), # Acetic Acid + Trimethylamine
        "Ethanol": ("Esterification", ["CC(=O)O", "CCO"]), # Acetic acid + Ethanol
        "Ethyl Acetate": ("Grignard_Reaction", ["CC(=O)OCC", "C[Mg]Br"]), # Ester + Grignard (rule needs adjustment for esters)
        "Toluene": ("Diels_Alder", ["CC1=CC=CC=C1", "C=C[CH]=O"]), # Toluene unlikely to undergo Diels-Alder as diene
        "Propane": ("Williamson_Ether_Synthesis", ["CCC[O-]", "CCBr"]) # Made up alkoxide for propane
    }

    for molecule_name, (reaction_name, reactants_smiles) in test_cases.items():
        print(f"Testing {molecule_name} with {reaction_name}:")
        products = apply_reaction(reaction_name, reactants_smiles)
        if products is not None:
            print(f"  Reactants: {reactants_smiles}")
            print(f"  Products: {products}")
        print("-" * 30)
