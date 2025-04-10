import utils
import iupac_to_smile
import preprocess_PC
import VSEPER
import inertie
import axe
import plan
import inversion
import subprocess
import sys

def extract_unique_bonds(graph):  #empeche d'avoir des liaisons appartir du graph redondante comme (1,2)=(2,1)
    seen = set()
    bonds = []
    for i_str, neighbors in graph.items():
        i = int(i_str)
        for j, _ in neighbors:
            key = tuple(sorted((i, j)))
            if key not in seen:
                seen.add(key)
                bonds.append(key)
    return bonds

def main():
    iupac_name = input("Entrez un nom IUPAC : ").strip()  #Entrée 
    smiles = iupac_to_smile.iupac_to_smiles(iupac_name)  #conversion par ""iupac_to_smile.py"""
    if not smiles:
        print(f"Erreur : Impossible de convertir '{iupac_name}' en SMILES.")
        return

    print(f"SMILES obtenu : {smiles}") #terminal debug
    molecule_data = preprocess_PC.process_molecule(smiles)#sauvegarde la sortie de preprocess_pc.py ( graph et matrice)
    molecule_data = {"iupac_name": iupac_name, **molecule_data}# ajout du IUPAC dans le fichier de save pour facilité la lecture

    if "masses" not in molecule_data:#debug de preprocess ?????
        print("Erreur : Les masses atomiques ne sont pas disponibles.")
        return

    coordinates = VSEPER.generate_3D_coordinates(smiles)# vseper.py pour les données 3D
    molecule_data["coordinates"] = coordinates.tolist()#sauvegarde de la sortie de vseper.py
    utils.save_data(molecule_data)

    inertia_axes = inertie.get_inertia_axes(molecule_data["coordinates"], molecule_data["masses"])
    molecule_data["inertia_axes"] = inertia_axes #calcul des axes d'inertie avec inertie.py a partir des données sauvegardé +sauvegarde
    utils.save_data(molecule_data)

    bonds = extract_unique_bonds(molecule_data["graph"])#traitement des donné du graphique moleculaire obtenus par preprocess.py
    detected_axes, tested_labels, candidate_axes = axe.detect_rotation_axes(  #axe.py pour testé les axes de symetrie
        molecule_data["coordinates"],
        molecule_data["masses"],
        bonds,
        molecule_data.get("inertia_axes", None)
    )
    molecule_data["rotation_axes"] = detected_axes#sauvegarde des axes de symetrie
    utils.save_data(molecule_data)

    mirror_planes = plan.detect_mirror_planes( #plan.py pour testé les plan et detecter une symetrie
        molecule_data["coordinates"],
        molecule_data["atom_symbols"],
        molecule_data["masses"],
        candidate_axes,
        logtest=False#DEBUG " true" pour voir quels plans ont été testé 
    )
    molecule_data["mirror_planes"] = mirror_planes#save
    utils.save_data(molecule_data)

    has_center = inversion.has_inversion_center(#inversion.py pour le centre d'inversion 
        molecule_data["coordinates"],
        molecule_data["atom_symbols"],
        molecule_data["masses"]
    )
    molecule_data["inversion_center"] = has_center#save
    utils.save_data(molecule_data)
#sortie affiché par le terminal


    print("\nAxes testés :")
    for label in tested_labels:
        print(f"  - {label}")

    print("\nAxes de rotation détectés :")
    for order, axes in detected_axes.items():
        print(f"  {order} :")
        for axis in axes:
            print(f"    - {axis['axis_label']}")

    print("\nPlans de symétrie détectés :")
    for plane in mirror_planes:
        print(f"  - {plane['plane_label']}")

    print("\nCentre d'inversion détecté :", "Oui" if has_center else "Non")

    subprocess.run([sys.executable, "visuel.py"])#debug de vseper.py
    subprocess.run([sys.executable, "visuel_axe.py"])#debug de intertie.py

if __name__ == "__main__":
    main()