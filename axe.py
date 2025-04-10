import numpy as np
from reperes import center_molecule
import utils

def get_rotation_matrix(axis, angle):
    axis = axis / np.linalg.norm(axis)#normalisation
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    ux, uy, uz = axis
    #matrice de rodrigues pour simuler la rotation d'angle Nd'un objet selon un axe de coordoné ux uy uz
    return np.array([
        [cos_theta + ux**2 * (1 - cos_theta), ux * uy * (1 - cos_theta) - uz * sin_theta, ux * uz * (1 - cos_theta) + uy * sin_theta],
        [uy * ux * (1 - cos_theta) + uz * sin_theta, cos_theta + uy**2 * (1 - cos_theta), uy * uz * (1 - cos_theta) - ux * sin_theta],
        [uz * ux * (1 - cos_theta) - uy * sin_theta, uz * uy * (1 - cos_theta) + ux * sin_theta, cos_theta + uz**2 * (1 - cos_theta)]
    ])

def detect_rotation_axes(coords, masses, bonds, inertia_axes=None, max_order=6, tolerance=0.3):
    coords = np.array(coords)
    masses = np.array(masses)
    coords, center_of_mass = center_molecule(coords, masses)#origine du repere

    detected_axes = {}
    candidate_axes = []
    tested_labels = []

    # Axes passant par les atomes (Centre/Atome)##### parcourt les coordonées des atomes excluant l'atome centrale 
    for i, atom_pos in enumerate(coords):
        if np.linalg.norm(atom_pos) > 0:
            label = f"Atom-{i+1}"
            axis = atom_pos / np.linalg.norm(atom_pos)#nomme chaque axe 
            candidate_axes.append((axis, label))
            tested_labels.append(label)

    # Axes passant par le milieu des liaisons (Centre-> centre liaison)
    seen = set()
    for i, j in bonds:
        key = tuple(sorted((i, j)))
        if key in seen:#anti doublons 
            continue
        seen.add(key) 

        mid_point = (coords[i] + coords[j]) / 2  #localise le milieu de liaison
        if np.linalg.norm(mid_point) > 0:
            label = f"MidBond-{i+1}-{j+1}"
            axis = mid_point / np.linalg.norm(mid_point)#nomme
            candidate_axes.append((axis, label))
            tested_labels.append(label)

    # Axes d'inertie (+ et -)
    if inertia_axes:
        for label, axis in inertia_axes.items():
            axis_vector = np.array(axis)
            for sign in [1, -1]:
                label_full = f"{label} ({'+' if sign == 1 else '-'})"
                candidate_axes.append((sign * axis_vector, label_full))
                tested_labels.append(label_full)

    # Détection 
    for axis, label in candidate_axes:
        detected_order, match_ratio = is_rotation_symmetric(coords, axis, range(max_order, 1, -1), tolerance)
        if detected_order:
            key = f"C{detected_order}"
            if key not in detected_axes:
                detected_axes[key] = []
            detected_axes[key].append({
                "axis_label": label,
                "axis_vector": axis.tolist(),
                "match_ratio": match_ratio
            })

        merged_axes = merge_similar_axes(detected_axes, merge_tolerance=0.1)#supprime les doublons
    return merged_axes, tested_labels, candidate_axes

def is_rotation_symmetric(coords, axis, orders, tolerance=0.3, min_match_ratio=1):#tolerance correspond à la distance autorisé entre 2 atomes avant et apres rotation## min_match_ratio correspond aux pourcentages d'atomes coinscident avant et apres rotation (ici 100%)
    for order in orders:
        angle = 2 * np.pi / order#theta 
        rotation_matrix = get_rotation_matrix(axis, angle)
        rotated_coords = np.dot(coords, rotation_matrix.T)
        distances = np.linalg.norm(coords - rotated_coords[:, None], axis=2)#comparaison avec les coordonné initial
        match_count = sum(np.min(distances, axis=1) <= tolerance)
        match_ratio = match_count / len(coords)
        if match_ratio >= min_match_ratio:
            return order, match_ratio
    return None, 0
## fonction pour eviter les doublons
def merge_similar_axes(detected_axes, merge_tolerance=0.1):#merge_tolerance EST Distance angulaire max pour considérer que deux vecteurs sont identiques ( evite les probleme de deformation)
    merged_axes = {}

    for order, axes in detected_axes.items():
        merged = []
        for axis_data in axes:
            axis = np.array(axis_data["axis_vector"])
            label = strip_direction_suffix(axis_data["axis_label"])
            match_ratio = axis_data["match_ratio"]

            # Vérifie si un axe équivalent existe déjà
            if not any(axes_equivalents(axis, np.array(existing["axis_vector"]), merge_tolerance) for existing in merged):
                merged.append({
                    "axis_vector": axis.tolist(),
                    "axis_label": label,
                    "match_ratio": match_ratio
                })

        merged_axes[order] = merged

    return merged_axes

def axes_equivalents(a, b, tol):
    return np.linalg.norm(a - b) < tol or np.linalg.norm(a + b) < tol

def strip_direction_suffix(label):
    return label[:-4] if label.endswith(" (+)") or label.endswith(" (-)") else label
