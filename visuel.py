import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection


ATOM_COLORS = ["b", "r", "g", "y", "orange", "purple", "c", "m"]

def load_molecule_data(json_path):
    
    with open(json_path, 'r') as file:
        data = json.load(file)
    return data

def plot_molecule(data):
    #Affiche la molécule en 3D à partir des données JSON
    coordinates = np.array(data['coordinates'])
    atom_symbols = [data['atom_symbols'][str(i)] for i in range(len(data['atoms']))]  
    graph = data['graph']
    
    # Attribution d'une couleur 
    unique_atoms = list(set(atom_symbols))
    atom_color_map = {atom: ATOM_COLORS[i % len(ATOM_COLORS)] for i, atom in enumerate(unique_atoms)}
    
    # Création de la figure
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Ajout des atomes
    for idx, (x, y, z) in enumerate(coordinates):
        color = atom_color_map[atom_symbols[idx]]
        ax.scatter(x, y, z, c=color, s=200, label=f"{atom_symbols[idx]}" if idx == 0 else "")

    # Ajout des liaisons
    edges = []
    for i, connections in graph.items():
        for j, _ in connections:
            edges.append([coordinates[int(i)], coordinates[int(j)]])
    ax.add_collection3d(Line3DCollection(edges, colors='k', linewidths=2))
    
    ###################
    ax.set_xlabel('X') #test
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Molécule : {data["iupac_name"]}')
    
    plt.legend()
    plt.show()

if __name__ == "__main__":
    json_path = "molecule_data.json"
    molecule_data = load_molecule_data(json_path)
    plot_molecule(molecule_data)
