#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;

typedef vector<vector<int>> Matrix;

void display_vector(double *vector, int L);
void display_matrix(double **matrix, int L, int C);
void display(Matrix lvl);
////////// create and delete

// fonction pour la creation d'un vecteur
double *newVector(int L);
int *newVectorInt(int L);
////////// on vector
// fonction qui sert a trouver un elemeent dans un vecteur

bool include(vector<int> const &list, int item);

class Cuthill
{
public:
    // Constructeur
    Cuthill(string filename);
    // Obtenir les voisins d'un nœud dans le graphe
    vector<int> getVoisin(int node);

    // Obtenir les voisins d'un nœud, en excluant les nœuds présents dans le vecteur "except"
    vector<int> getVoisin(int node, vector<int> except);

    // Générer la structure de niveau du nœud
    Matrix generateStruct(int node);

    // Appliquer la factorisation sur le profil optimisé
    void factorize();

    // Utilisé lorsque nous résolvons le système factorisé
    void resolutionInf();
    void resolutionSup();

    // Résoudre le système en utilisant l'algorithme de Cuthill-McKee
    void solve();

    // Afficher le résultat
    void displayResult();

    // Calculer l'erreur moyenne
    double meanError();
    // Setter pour le profil optimisé de A
    void ChangeAp(int i, int j, double value);
    // Trouver le premier sommet pour commencer l'algorithme de Cuthill-McKee
    int getFirstNoeud();

    // Utiliser l'algorithme de Cuthill-McKee et stocker le renumérotage des sommets dans le nœud
    void CuthillMckee(int first_node);

    // Optimiser le renumérotage des sommets avec l'algorithme de Cuthill-McKee inverse
    void cuthillMckceInverse(int first_node);

    // Transformer notre profil en profil optimisé (AOptimise)
    void optimizeProfil();
    // Getter pour le profil de A (avant l'optimisation)
    double getAP(int i, int j);

    // Getter pour le profil optimisé de A
    double FactorisationAP(int i, int j);

private:
    vector<double> A;         // Le profil de A
    vector<double> AOptimise; // Le profil optimisé de A (optimisé car dim(AOptimise) < dim(A))

    double *b;   // Second membre
    double *bp;  // Résultat de P.double * b (où P est la matrice de permutation du graphe de A au graphe de AOptimise)
    double *x;   // Solution du système
    double *xp;  // Résultat de AOptimise * xp = bp
    double *sol; // La vraie solution

    int *nDiag;         // nDiag (1) associé au profil A
    int *nDiagOptimise; // nDiag (1) associé au profil AOptimise
    int *p;             // Liste du premier terme non nul sur chaque colonne de A
    int *pOptimise;     // Liste du premier terme non nul sur chaque colonne de AOptimise
    int *nodes;         // Contient la permutation de chaque nœud (ancien nœud en tant qu'indice -> nouveau nœud en tant que valeur)
    int size;           // La taille de la matrice
};

Cuthill::Cuthill(string filename)
{
    string temp("");
    double value(0);
    int profil_size(0);
    bool test(false);
    ifstream stream(filename);

    if (stream)
    {
        stream >> size;
        // init data
        nDiag = newVectorInt(size);
        nDiagOptimise = newVectorInt(size);
        p = newVectorInt(size);
        pOptimise = newVectorInt(size);

        nodes = newVectorInt(size);
        b = newVector(size);
        bp = newVector(size);
        x = newVector(size);
        xp = newVector(size);
        sol = newVector(size);

        // store profil in A
        for (int i = 0; i < size; i++)
        {
            test = false;

            for (int j = 0; j < size; j++)
            {
                stream >> value;
                if (j <= i)
                {
                    if (test == false && value != 0)
                    {
                        test = true;
                    }
                    if (test == true)
                    {
                        A.push_back(value);
                        profil_size++;
                    }
                }
            }
            // compute nDiag and p
            nDiag[i] = profil_size - 1;
            p[i] = i - ((i > 0) ? (nDiag[i] - nDiag[i - 1] - 1) : 0);
        }
        cout << endl
             << "++++++++++++++++++++++++ Matrice A +++++++++++++" << endl
             << endl;
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                if (j < p[i])
                {
                    cout << "* ";
                }
                else
                {
                    cout << getAP(i, j) << " ";
                }
            }
            cout << endl;
        }

        // read second member
        for (int i = 0; i < size; i++)
        {
            stream >> b[i];
            nodes[i] = -1;
        }
        cout << endl
             << "Second member   b : ";
        display_vector(b, size);
    }
    else
    {
        cout << "Le fichier n'existe pas" << endl;
        exit(1);
    }
}

// get voisins of node in graph
// we get it in profil of matrix
vector<int> Cuthill::getVoisin(int node)
{

    vector<int> voisin;

    for (int i = 0; i < size; i++)
    {
        if (i == node)
            continue;
        else if (i < node)
        {
            if (getAP(node, i) != 0)
                voisin.push_back(i);
        }
        else
        {
            if (getAP(i, node) != 0)
                voisin.push_back(i);
        }
    }

    return voisin;
}
// prendre les voisin de noeud a part les neoud dans except
vector<int> Cuthill::getVoisin(int node, vector<int> except)
{

    vector<int> voisin;

    for (int i = 0; i < size; i++)
    {
        if (i == node)
            continue;
        else if (i < node)
        {
            if (getAP(node, i) != 0 && !include(except, i))
                voisin.push_back(i);
        }
        else
        {
            if (getAP(i, node) != 0 && !include(except, i))
                voisin.push_back(i);
        }
    }

    return voisin;
}
/**
 * Génère la structure de niveau du nœud.
 * La structure de niveau du nœud est une liste des voisins du nœud en suivant la distance entre eux.
 * Ainsi, nous pouvons trouver le nœud le plus éloigné de ce nœud dans cette structure.
 */
Matrix Cuthill::generateStruct(int node)
{

    Matrix lvl;
    vector<int> selected;
    int k = 0;

    lvl.push_back(vector<int>(1, node));
    selected.push_back(node);

    while (selected.size() != size)
    { // while we haven't store all node in lvl
        vector<int> voisins;
        for (int current_node : lvl[k])
        { // for each voisin
            vector<int> voisin_of_voisins = getVoisin(current_node, selected);
            for (int voisin : voisin_of_voisins)
            {
                // we store all voisin in selected
                selected.push_back(voisin);
                voisins.push_back(voisin);
            }
        }
        lvl.push_back(voisins);
        k++;
    }

    return lvl;
}

int Cuthill::getFirstNoeud()
{

    Matrix lvl_first_node, lvl_current_node;
    vector<int> nodeAlreadyComputed;
    int first_node = 0;                       // (1) in lesson
    int e_first_node = 0, e_current_node = 0; // eccentricities
    bool first_node_found = false;

    while (!first_node_found)
    {
        lvl_first_node = generateStruct(first_node); // (2) in lesson
        e_first_node = lvl_first_node.size() - 1;
        first_node_found = true; // I assume that we found the node here

        for (int current_node : lvl_first_node[lvl_first_node.size() - 1])
        { // (3) in lesson
            if (include(nodeAlreadyComputed, current_node))
            {
                continue;
            }
            lvl_current_node = generateStruct(current_node);
            e_current_node = lvl_current_node.size() - 1;
            // if we find node that is the most eccentric that n.
            if (e_current_node > e_first_node)
            {
                // we save that node
                nodeAlreadyComputed.push_back(first_node);
                first_node = current_node;
                first_node_found = false;
                break;
            }
            else
            {
                nodeAlreadyComputed.push_back(current_node);
            }
        }
    }
    return first_node;
}

void Cuthill::CuthillMckee(int first_node)
{

    vector<int> fixed;

    int current_index = 0,
        node_queue[size];

    node_queue[0] = first_node;
    nodes[node_queue[0]] = current_index;
    fixed.push_back(node_queue[0]);

    for (int cn = 0; cn < size; cn++)
    {
        vector<vector<int>> dict;
        vector<int> voisins = getVoisin(node_queue[cn], fixed);

        if (voisins.size() > 0)
        { // si on trouve le voisin

            for (int voisin : voisins)
            { // stocker le voisin et son nombre de voisin
                dict.push_back({(int)getVoisin(voisin, fixed).size(),
                                voisin});
            }
            sort(dict.begin(), dict.end(), [](vector<int> const &a, vector<int> const &b)
                 { return a[0] < b[0]; });
            for (int i = 0; i < dict.size(); i++)
            {
                node_queue[current_index + (i + 1)] = dict[i][1];
                nodes[dict[i][1]] = current_index + (i + 1);
                fixed.push_back(dict[i][1]);
            }
            current_index += dict.size();
        }
    }
}

void Cuthill::cuthillMckceInverse(int first_node)
{
    CuthillMckee(first_node);
    for (int i = 0; i < size; i++)
        nodes[i] = (size - 1) - nodes[i];
}

void Cuthill::optimizeProfil()
{

    double value = 0;
    int profil_size = 0;
    bool test = false;
    for (int i = 0; i < size; i++)
    {
        test = false;
        for (int j = 0; j <= i; j++)
        {
            value = getAP(nodes[i], nodes[j]);
            if (test == false && value != 0)
                test = true;

            if (test == true)
            {
                AOptimise.push_back(value);
                profil_size++;
            }
        }
        // we initialize nDiag and p of AOptimise
        nDiagOptimise[i] = profil_size - 1;
        pOptimise[i] = i - ((i > 0) ? (nDiagOptimise[i] - nDiagOptimise[i - 1] - 1)
                                    : 0);
    }

    cout << endl
         << "++++++++++++++++++++++++ Matrice réarrangée +++++++++++++" << endl;
    cout << endl
         << "++++++++++++++++++++++++ Profil A optimisé +++++++++++++" << endl
         << endl;
    ofstream file("cuthil.txt");
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            if (j < pOptimise[i])
            {
                cout << ". ";
                file << ". ";
            }
            else
            {
                cout << FactorisationAP(i, j) << " ";
                file << FactorisationAP(i, j) << " ";
            }
        }
        file << endl;
        cout << endl;
    };
    cout << endl
         << "Profil A est stocké dans un fichier cuthil.txt" << endl;
    double sum = 0;
    // we compute b' here
    for (int i = 0; i < size; i++)
    {
        sum = 0;
        for (int k = 0; k < size; k++)
        {
            if (k == nodes[i])
                sum += b[k];
        }
        bp[i] = sum;
    }
}

/**

Facatorisation en LDLT
 */

void Cuthill::factorize()
{

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < i; j++)
        {
            for (int k = 0; k <= j - 1; k++)
                ChangeAp(i, j, FactorisationAP(i, j) - (FactorisationAP(i, k) * FactorisationAP(k, k) * FactorisationAP(j, k)));
            ChangeAp(i, j, FactorisationAP(i, j) / FactorisationAP(j, j));
        }

        for (int k = 0; k <= i - 1; k++)
            ChangeAp(i, i, FactorisationAP(i, i) - FactorisationAP(k, k) * FactorisationAP(i, k) * FactorisationAP(i, k));
    }
}

void Cuthill::resolutionInf()
{
    double sum(0);
    // 1- compute x where L * x = b
    for (int k = 0; k < size; k++)
    {
        if (k == 0)
        {
            xp[k] = bp[k];
        }
        else
        {
            sum = 0;
            for (int j = 0; j <= k; j++)
                sum += FactorisationAP(k, j) * xp[j];
            xp[k] = bp[k] - sum;
        }
    }

    for (int i = 0; i < size; i++)
    {
        xp[i] /= FactorisationAP(i, i);
    }
}

void Cuthill::resolutionSup()
{
    double sum(0);
    for (int k = size - 1; k >= 0; k--)
    {
        if (k != size - 1)
        {
            sum = 0;
            for (int j = k + 1; j < size; j++)
                sum += FactorisationAP(j, k) * xp[j];
            xp[k] = xp[k] - sum;
        }
    }
}

/**
 * Resolution X = P*X'
 */
void Cuthill::solve()
{

    double sum(0);

    for (int i = 0; i < size; i++)
    {
        sum = 0;
        for (int k = 0; k < size; k++)
        {
            if (i == nodes[k])
                sum += xp[k];
        }
        x[i] = sum;
    }
}

void Cuthill::displayResult()
{
    const double error = meanError();
    cout << endl
         << "++++++++++++++ Permurtation des noeud +++++++++++++++++===" << endl
         << endl;
    for (int i = 0; i < size; i++)
        cout << " " << i;
    cout << endl;
    for (int i = 0; i < size; i++)
        cout << " " << ((nodes[i] == -1) ? "." : to_string(nodes[i]));
    cout << endl
         << endl
         << "Solution trouvé: ";
    display_vector(x, size);
    cout << endl
         << endl;
    cout << "Erreur entre la solution exate et la solution trouvé : " << error << endl;
    cout << endl;
}

/**
Calculer l'erreur entre la solution exacte et la solution à calculer
 */
double Cuthill::meanError()
{
    double error(0);
    for (int i = 0; i < size; i++)
    {
        error += abs(x[i] - sol[i]);
    }
    error /= size;
    return error;
}

double Cuthill::getAP(int i, int j)
{
    if (j >= p[i] && j <= i)
        return A[nDiag[i] - i + j];
    else if (i >= p[j] && i < j)
        return A[nDiag[j] - j + i];
    else
        return 0;
}

double Cuthill::FactorisationAP(int i, int j)
{
    if (j >= pOptimise[i] && j <= i)
        return AOptimise[nDiagOptimise[i] - i + j];
    else
        return 0;
}

void Cuthill::ChangeAp(int i, int j, double value)
{
    if (j >= pOptimise[i])
        AOptimise[nDiagOptimise[i] - i + j] = value;
}

/***************************Genere les affichge***************************/

void display_vector(double *vector, int L)
{
    for (int i = 0; i < L; i++)
    {
        cout
            << "  "
            << setw(2) << setprecision(15)
            << vector[i] << " ";
    }
    cout << endl;
}

void display_matrix(double **matrix, int L, int C)
{
    cout << "\n";
    for (int i = 0; i < (int)L; i++)
    {
        cout << "[";
        for (int j = 0; j < (int)C; j++)
        {
            cout
                << "  "
                << setw(2) << setprecision(15)
                << matrix[i][j] << ",";
        }
        cout << "  ],\n";
    }
}

void display(Matrix lvl)
{

    cout << "{";
    for (int i = 0; i < lvl.size(); i++)
    {
        cout << "\n   N" << i << " : { ";
        for (int j = 0; j < lvl[i].size(); j++)
        {
            cout << lvl[i][j] + 1 << " ";
        }
        cout << "},";
    }
    cout << "\n}" << endl;
}

double *newVector(int L)
{
    double *v;
    v = new double[L];
    if (v == nullptr)
        exit(1);
    return v;
}
int *newVectorInt(int L)
{
    int *v;
    v = new int[L];
    if (v == nullptr)
        exit(1);
    return v;
}

bool include(vector<int> const &list, int item)
{
    return find(list.begin(), list.end(), item) != list.end();
}