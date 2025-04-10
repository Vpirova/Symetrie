import numpy as np
from reperes import center_molecule

def reflect(coords, normal):
    normal = normal / np.linalg.norm(normal)
    projection = np.dot(coords, normal)[:, np.newaxis]  # (N, 1)
    return coords - 2 * projection * normal  # (N, 3)

def planes_equivalents(n1, n2, tol):
    return np.linalg.norm(n1 - n2) < tol or np.linalg.norm(n1 + n2) < tol

def merge_similar_planes(detected_planes, merge_tolerance=0.1):
    merged = []
    for plane in detected_planes:
        n = np.array(plane["normal"])
        label = plane["plane_label"]

        if not any(planes_equivalents(n, np.array(p["normal"]), merge_tolerance) for p in merged):
            merged.append({
                "plane_label": label,
                "normal": n.tolist()
            })
    return merged

def generate_planes_from_two_vectors(candidate_axes):
    """
    Construit des plans définis par deux vecteurs testés.
    La normale du plan est donnée par le produit vectoriel v1 × v2.
    Le plan résultant est celui qui CONTIENT les vecteurs v1 et v2.
    """
    new_planes = []
    for i in range(len(candidate_axes)):
        for j in range(i + 1, len(candidate_axes)):
            v1, label1 = candidate_axes[i]
            v2, label2 = candidate_axes[j]
            v1 = np.array(v1)
            v2 = np.array(v2)
            normal = np.cross(v1, v2)
            if np.linalg.norm(normal) < 1e-6:
                continue  # vecteurs presque colinéaires → pas de plan défini
            label = f"Plane({label1}, {label2})"
            new_planes.append((normal, label))
    return new_planes

def detect_mirror_planes(coords, atom_symbols, masses, candidate_axes, tolerance=0.2, logtest=False):
    mirror_planes = []

    coords, _ = center_molecule(np.array(coords), np.array(masses))

    all_planes = generate_planes_from_two_vectors(candidate_axes)

    for normal_vector, label in all_planes:
        normal = np.array(normal_vector) / np.linalg.norm(normal_vector)
        reflected_coords = reflect(coords, normal)

        unmatched = set(range(len(coords)))
        match_count = 0
        fixed_count = 0

        while unmatched:
            i = unmatched.pop()
            ri = reflected_coords[i]
            sym_i = atom_symbols[str(i)]
            found = False
            for j in list(unmatched):
                if atom_symbols[str(j)] != sym_i:
                    continue
                if np.linalg.norm(ri - coords[j]) <= tolerance:
                    unmatched.remove(j)
                    match_count += 1
                    found = True
                    break
            if not found:
                # test si l'atome est invariant par réflexion (sur le plan)
                if np.linalg.norm(ri - coords[i]) <= tolerance:
                    fixed_count += 1
                else:
                    break  # un atome non apparié ni fixé → on rejette ce plan

        total = len(coords)
        matched_atoms = match_count * 2 + fixed_count

        if logtest:
            print(f"Test du plan {label} : {matched_atoms}/{total} atomes traités")

        if matched_atoms == total:
            mirror_planes.append({
                "plane_label": label,
                "normal": normal.tolist()
            })

    return merge_similar_planes(mirror_planes)