# Code Julia pour le projet REOP

## Utilisation

Ce dossier contient un code Julia permettant de bien démarrer le projet du cours de recherche opérationnelle. Il peut bien sûr être utilisé en Julia, mais aussi appelé via Python grâce à la librarie [pyjulia](https://pyjulia.readthedocs.io/en/stable/). Si vous débutez en Julia, n'hésitez pas à vous référer à mon tutoriel https://github.com/gdalle/IntroJulia.

Lors de la première exécution, il est probable que vous ne possédiez pas tous les packages Julia nécessaires. Dans ce cas, il suffit de suivre les instructions données par la console et d'installer chacun graĉe aux commandes suivantes :

```julia
import Pkg
Pkg.add("packagename")
```

S'il vous manque des packages Python (notamment folium), c'est dans le fichier `plot.jl` qu'il faut les installer via Conda, en décommentant la ligne appropriée au début du fichier.

Le notebook Jupyter `Projet REOP.ipynb` présente un exemple d'utilisation des principales fonctions.

## Contenu

Toutes les fonctions sont importées par le fichier `import_all.jl`. Voici une brève description des autres fichiers :

- `cost.jl` : Calcul du coût d'une solution.
- `emballage.jl` : Définition de la classe `Emballage` et lecture à partir d'une chaîne de caractères.
- `feasibility.jl` : Satisfaction des contraintes par une solution.
- `fournisseur.jl` : Définitin de la classe `Fournisseur` et lecture à partir d'une chaîne de caractères.
- `graphe.jl` : Définition de la classe `Graphe` (qui stocke les distances entre sites) et lecture à partir d'une chaîne de caractères.
- `instance.jl` : Définition de la classe `Instance`, qui regroupe tous les paramètres d'une instance ainsi qu'une solution (vide par défaut), et lecture à partir d'un fichier.
- `plot.jl` : Outils de visualisation d'une instance ou d'une solution.
- `route.jl` : Définition de la classe `Route` et lecture à partir d'une chaîne de caractères.
- `solution.jl` : Lecture d'une solution à partir d'un fichier, et stockage à l'intérieur d'une instance.
- `usine.jl` : Définitin de la classe `Usine` et lecture à partir d'une chaîne de caractères.
- `write.jl` : Écriture des instances et solutions dans des fichiers texte.
