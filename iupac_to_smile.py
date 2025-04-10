import ssl
from pubchempy import get_compounds
# Désactiver la vérification SSL 
ssl._create_default_https_context = ssl._create_unverified_context #test

def iupac_to_smiles(iupac_name):
   
    if not iupac_name:
        raise ValueError("Le nom IUPAC est vide.")

    try:
        compound = get_compounds(iupac_name, 'name')#interroge la base de donné de pubchem pour passer de iupac à smile 
        if compound:
            return compound[0].isomeric_smiles
        else:
            return None
    except Exception as e:
        print(f" Erreur lors de la conversion de {iupac_name} : {e}")
        return None

if __name__ == "__main__":
    iupac_name = input("Entrez un nom IUPAC, sans majuscule: ").strip()
    smiles = iupac_to_smiles(iupac_name)
    if smiles:
        print(f"{iupac_name} -> {smiles}")# sortie terminal
    else:
        print(" Conversion échouée.")
