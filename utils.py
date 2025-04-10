import json
#gestion du json

#création du fichier
def save_data(data, filename="molecule_data.json"):
    
    with open(filename, "w") as f:
        json.dump(data, f, indent=4)

#édition du fichier
def load_data(filename="molecule_data.json"):
    
    try:
        with open(filename, "r") as f:
            return json.load(f)
    except FileNotFoundError:
        print(f" Fichier {filename} introuvable. Retour d'un objet vide.")
        return {}  