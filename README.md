# Algorithme Cuthill Mackee


## Soit le systeme (S)

## objectif

##### On veut résoudre un système d'équations linéaires (S) : Ax = b sachant que A

##### est creuse.

##### Pour cela, on se propose de stocker la matrice par la méthode de rangement

##### profil en utilisant un profil

##### optimisé grâce à la renumérotation de sommets pour ensuite résoudre le

##### système avec la méthode de

##### factorisation LDL_t.

##### Pour la renumérotation, on va utiliser l'algorithme de Cuthill McKee inverse.

## Algorithme

#### Pour optimiser le profil, on renumérote les sommets du graphe comme

#### suit:

#### ◦ choix du premier noeud avec l'algorithme de Cuthill McKee

#### ◦ les autres sommets sont renumérotés successivement en comptant les

#### voisins non renumérotés du

#### dernier sommet renuméroté

#### Algorithme de Cuthill McKee inverse

#### L'algorithme de Cuthill McKee permet de choisir de façon plus ou moins

#### optimale premier noeud à

#### renuméroter.

#### L'algorithme consiste en 4 étapes:

#### 1.choix du noeud à renuméroter n: on a choisi le noeud 0;


#### 2.Tant qu'on n'a pas épuisé tous les sommets, on répète les étapes 3. et 4. :

#### 3.calculer e(n);

#### 4.choix de s dans N e(n)

#### ▪ calcul de e(s)

#### ▪ si e(s) > e(n), on revient à 3.

#### ▪ sinon, on revient à 4. en passant à l'élément suivant

#### Une fois le premier noeud choisi, on renumérote le reste des sommets

#### comme suit:

#### Pour chaque sommet renuméroté, on cherche ses voisins et on attribue le

#### numéro suivant selon les voisins

#### ayant le moins de voisins non renumérotés.

#### L'algorithme de Cuthill McKee inverse consiste ensuite à appliquer la

#### formule suivante pour obtenir la

#### renumérotation finale:

#### Si i est le sommet à renuméroter ey σ(i) sa transformée par l'algorithme de

#### Cuthill McKee, on a:

#### σ ′ (i) = DIM − 1 − σ(i) avec DIM: la dimension de la matrice

#### exemple

#### Prenons le sommet s = 9 et cherchons son exentricité.

#### Rangement profil

#### Après la permutation, on stocke les éléments du profil de A' dans le tableau

#### AP et les indices, dans AP des


#### éléments de la diagonale dans le tableau nDiag.

#### La répartition des éléments de A dans les tableaux AP et nDiag est la

#### suivante:

```
AP [ nDiag [ i ]− i + j ]=
```
### ( AP [ nDiag [ i ]− i + j ]−∑ k = i

```
j
( AP [ nDiag [ i ]− i + k ]∗ AP [ nDiag [ k ]]∗ AP [ nDiag [ j ]− j + k ]))
( AP [ nDiag [ j ]])
```
### AP [ nDiag [ i ]]= AP [ nDiag [ i ]]− k ∑= p

```
i
```
```
i
AP [ nDiag [ k ]]∗( AP [ nDiag [ i ]− i + k ])^2
```
#### Résolution de (s)

#### La résolution de (S) s'e ectue en appliquant la factorisation LDL_t:ff

#### A = L. D. Lt avec:

#### L[i][j] = AP [nDiag[i] − i + j] ; j < i

#### D[i] = AP [nDiag[i]]


