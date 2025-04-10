from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
#### coordonné 3D
def generate_3D_coordinates(smiles, num_confs=100):
    #Génère une conformation 3D 
    if not smiles:
        raise ValueError(" Le SMILES est vide ou invalide.")
    
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if mol is None:
        raise ValueError(f" Impossible de générer la molécule à partir de {smiles}.")
    
    params = AllChem.ETKDGv3()
    params.randomSeed = 42# seed fixé pour rendre le debug plus simple
    
    # Générer plusieurs conformations et prendre la plus stable(100) pour evité des cycles trop déformé entre bateau et chaise 
    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
    if not conf_ids:
        raise ValueError(" Échec de l'intégration 3D des conformères.")
    
    best_conf_id = min(conf_ids, key=lambda conf_id: AllChem.UFFGetMoleculeForceField(mol, confId=conf_id).CalcEnergy())# utilise UFF pour calculer l'energie des 100 conformere et prend la plus stable
    ff_uff = AllChem.UFFGetMoleculeForceField(mol, confId=best_conf_id)#ajustement ( raffinement energetique ) , meilleurs resultats ? 
    ff_uff.Minimize()#OPTIMIsation local, meilleurs resultats ?? 
    
    conf = mol.GetConformer(best_conf_id)#sortie de la meilleure version 
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])#sauvegarde des coordonnées d'atome
    
    # Correction des artefacts ???? a faire ? 
    #coords = correct_coordinate_artifacts(coords)
    
    return coords



if __name__ == "__main__":#debug
    smiles = input("Entrez un SMILES pour tester la génération 3D : ").strip()
    try:
        coords = generate_3D_coordinates(smiles)
        print(" Coordonnées 3D corrigées :\n", coords)
    except ValueError as e:
        print(f" Erreur : {e}")
