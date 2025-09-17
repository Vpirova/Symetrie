Projet Python – SMILE 
PIROVANI
Date : [10/04/2025]

Description :
-------------
Ce programme permet d’analyser les éléments de symétrie d’une molécule à partir de son nom IUPAC. Il convertit le nom en SMILES, génère une géométrie 3D approximative du conformère le plus stable, puis identifie les axes de rotation, plans de symétrie et centres d'inversion.

Utilisation :
-------------
1. Lancer `main.py`
2. Entrer un nom IUPAC valide (ex. "methane", "benzene", etc.)
3. Observer les résultats affichés dans le terminal et les visualisations générées

Fichiers du projet :
--------------------
- `main.py` : script principal à exécuter
- `iupac_to_smile.py` : conversion IUPAC → SMILES
- `preprocess_PC.py` : génération du graphe moléculaire et des données de base
- `VSEPER.py` : estimation des coordonnées 3D
- `inertie.py` : calcul des axes d’inertie
- `axe.py` : détection des axes de rotation
- `plan.py` : détection des plans de symétrie
- `inversion.py` : test du centre d’inversion
- `utils.py` : sauvegarde/chargement des données
- `visuel.py`, `visuel_axe.py` : scripts de visualisation

Dépendances :
-------------
Voir `biblio.txt` – installation rapide avec :
    pip install -r biblio.txt

    

N.B 17/09/2025 : 
-------------
![status](https://img.shields.io/badge/Status-PANNE%20CONNUE-red)

> 🔴 **Origine — Conversion IUPAC → SMILES**  
> La conversion repose sur PubChemPy (requêtes en ligne vers PubChem). En cas de problème réseau/SSL ou de réponse vide coté serveur, la fonction renvoie `Erreur : Impossible de convertir '' en SMILES.` et le script s’arrête.  Un correctif  ( petit dictionnaire local) sera intégré plus tard.

