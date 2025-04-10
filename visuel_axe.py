import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from reperes import center_molecule
from inertie import get_inertia_axes

def load_molecule_data(json_path):
    
    with open(json_path, 'r') as file:
        data = json.load(file)
    return data

def plot_molecule_with_axes(data, scale=2):
    #affichage
    coordinates = np.array(data['coordinates'])
    masses = np.array(data['masses'])
    graph = data['graph']
    
    
    coordinates, center_of_mass = center_molecule(coordinates, masses)

   
    inertia_axes = get_inertia_axes(coordinates, masses)

    # Création de la figure 3D
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Affichage des atomes
    ax.scatter(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2], c='b', s=100, label="Atomes")
    
    # Affichage des liaisons
    edges = []
    for i, connections in graph.items():
        for j, _ in connections:
            edges.append([coordinates[int(i)], coordinates[int(j)]])
    ax.add_collection3d(Line3DCollection(edges, colors='k', linewidths=2))
    
    # **Affichage des axes d'inertie centrés sur le centre de masse**
    colors = ["r", "g", "b"]
    labels = ["Principal Axis Z", "Principal Axis X", "Principal Axis Y"]
    
    for idx, axis in enumerate(inertia_axes.values()):
        axis = np.array(axis) * scale  # Mise à l'échelle pour visibilité
        ax.plot(
            [center_of_mass[0], center_of_mass[0] + axis[0]],
            [center_of_mass[1], center_of_mass[1] + axis[1]],
            [center_of_mass[2], center_of_mass[2] + axis[2]],
            linestyle="--",
            linewidth=2,
            color=colors[idx],
            label=labels[idx]
        )

   
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f"Molécule : {data['iupac_name']} et Axes d'Inertie")
    ax.legend()
    plt.show()

if __name__ == "__main__":
    json_path = "molecule_data.json"
    molecule_data = load_molecule_data(json_path)
    plot_molecule_with_axes(molecule_data)
