import numpy as np

def compute_inertia_tensor(coords, masses):
    #tenseur
    if len(coords) != len(masses):
        raise ValueError("Nombre d'atomes et de masses incohérent.")

    I = np.zeros((3, 3))
    
    for i in range(len(coords)):
        x, y, z = coords[i]
        m = masses[i]
        
        # Éléments diagonaux
        I[0, 0] += m * (y**2 + z**2)
        I[1, 1] += m * (x**2 + z**2)
        I[2, 2] += m * (x**2 + y**2)

        # Éléments hors-diagonaux
        I[0, 1] -= m * x * y
        I[0, 2] -= m * x * z
        I[1, 2] -= m * y * z

    I[1, 0] = I[0, 1]
    I[2, 0] = I[0, 2]
    I[2, 1] = I[1, 2]

    return I

def get_principal_axes(coords, masses):
    
    I = compute_inertia_tensor(coords, masses)#on part du tenseur I
    
    eigenvalues, eigenvectors = np.linalg.eigh(I)  # on prend la matrice I et on en extrait ses valeurs propres.
    sorted_indices = np.argsort(eigenvalues)[::-1]  # Tri pour debug z , y ,x 
    sorted_axes = eigenvectors[:, sorted_indices]

    return sorted_axes


def get_inertia_axes(coords, masses):#convertis les  vecteurs propres contenus dans une matrice rn listes
    
    principal_axes = get_principal_axes(coords, masses)

    return {
        "principal_axis_Z": principal_axes[:, 0].tolist(),
        "principal_axis_X": principal_axes[:, 1].tolist(),
        "principal_axis_Y": principal_axes[:, 2].tolist()
    }

if __name__ == "__main__":
    import utils  #  test
    data = utils.load_data()
    
    if "coordinates" in data and "masses" in data:
        inertia_axes = get_inertia_axes(np.array(data["coordinates"]), np.array(data["masses"]))
        print("Axes d'inertie calculés :", inertia_axes)
    else:
        print(" Erreur : Données manquantes dans `molecule_data.json`")
