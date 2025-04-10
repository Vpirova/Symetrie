import numpy as np
from reperes import center_molecule

def has_inversion_center(coords, atom_symbols, masses, tolerance=0.2):
    coords = np.array(coords)
    masses = np.array(masses)
    coords, _ = center_molecule(coords, masses)  # Centrage sur le centre de masse

    unmatched = set(range(len(coords)))  # Atomes non encore appariés

    while unmatched:
        i = unmatched.pop()
        ri = coords[i]
        symbol_i = atom_symbols[str(i)]
        found_pair = False

        for j in list(unmatched):
            rj = coords[j]
            symbol_j = atom_symbols[str(j)]
            if symbol_i != symbol_j:
                continue
            if np.linalg.norm(ri + rj) <= tolerance:
                unmatched.remove(j)
                found_pair = True
                break

        if not found_pair:
            return False  # Pas de paire opposée trouvée

    return True
