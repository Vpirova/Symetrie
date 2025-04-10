import numpy as np

def get_atomic_masses(atom_symbols):
    #ableau des masses atomiques + symboles chimiques
    atomic_masses = {
        "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011,
        "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305,
        "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06, "Cl": 35.45, "Ar": 39.948,
        "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996,
        "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38,
        "Ga": 69.723, "Ge": 72.63, "As": 74.922, "Se": 78.96, "Br": 79.904, "Kr": 83.798,
        "Rb": 85.468, "Sr": 87.62
    }
    
    return np.array([atomic_masses.get(atom_symbols[str(i)], 12.011) for i in range(len(atom_symbols))]) #associe les symboles présent dans la molécules avec les valeurs du tableau 

def center_molecule(coords, masses):
    #Centre la molécule sur son centre de masse
    total_mass = np.sum(masses)
    #moyenne pondéré 
    center_of_mass = np.sum(coords * masses[:, np.newaxis], axis=0) / total_mass
    return coords - center_of_mass, center_of_mass

if __name__ == "__main__":
    import utils  
    data = utils.load_data()
    
    if "coordinates" in data and "atom_symbols" in data:
        masses = get_atomic_masses(data["atom_symbols"])
        new_coords, center = center_molecule(np.array(data["coordinates"]), masses)
        
        print(f" Centre de masse calculé : {center}")
        print(f" Nouvelles coordonnées centrées :\n{new_coords}")
    else:
        print(" Erreur : Données manquantes dans `molecule_data.json`")
