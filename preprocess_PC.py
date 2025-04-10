from collections import defaultdict
import numpy as np
from rdkit import Chem
import reperes



def load_molecule_from_smiles(smiles):
    
    if not smiles:
        raise ValueError("Le SMILES est vide ou invalide.")

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles)) #rdkit crée la molécule à partir de smile avec hydrogene
    if mol is None:
        raise ValueError(f" molecule non valide pour le SMILES : {smiles}")
    
    return mol

def create_atoms(mol):  #pour avoir le symbole associé a chaque atome
    #liste
    atom_dict = {}
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    
  
    for a in mol.GetAromaticAtoms():#pointe les C aromatique
        i = a.GetIdx()
        atoms[i] = (atoms[i], 'aromatic')
    
    
    atoms = [atom_dict.setdefault(a, len(atom_dict)) for a in atoms]
    
    return np.array(atoms)

def create_ijbonddict(mol):
    #graph moléculaire
    bond_dict = {}
    i_jbond_dict = defaultdict(list)
    
    for b in mol.GetBonds():#parcours les liaisons
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()#associé les 2 atomes d'une mm liaisons
        bond_type = str(b.GetBondType())
        bond = bond_dict.setdefault(bond_type, len(bond_dict))#single double triple 
        
        i_jbond_dict[i].append((j, bond))#le dico des liaisons 
        i_jbond_dict[j].append((i, bond))
    
    return i_jbond_dict

def get_adjacency_matrix(mol):
    #matrice d'adjacence par rdkit
    return Chem.GetAdjacencyMatrix(mol)

def process_molecule(smiles):#fonction principale pour remplir le fichier de sauvegarde des données 
    
    mol = load_molecule_from_smiles(smiles)

    atoms = create_atoms(mol)
    i_jbond_dict = create_ijbonddict(mol)
    adjacency_matrix = get_adjacency_matrix(mol)

    atom_symbols = {str(i): a.GetSymbol() for i, a in enumerate(mol.GetAtoms())}  
    
    # Ajout du calcul des masses ################### dans reperes.py il y a le tableau périodique 
    masses = reperes.get_atomic_masses(atom_symbols)

    return {
        "smiles": smiles,
        "atoms": atoms.tolist(),
        "atom_symbols": atom_symbols,
        "masses": masses.tolist(),  
        "graph": {k: v for k, v in i_jbond_dict.items()},
        "adjacency_matrix": adjacency_matrix.tolist()
    }

if __name__ == '__main__':#debug
    smiles = input("Entrez un SMILES : ").strip()
    
    try:
        molecule_data = process_molecule(smiles)
        print(" Molécule traitée avec succès :", molecule_data)
    except ValueError as e:
        print(f" Erreur : {e}")
