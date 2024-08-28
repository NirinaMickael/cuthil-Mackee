#include <iostream>
#include <vector>   
#include <fstream>    
#include <cmath>      
#include <iomanip>    
#include <algorithm>  
#include"Cuthill.cpp"
using namespace std;
int main() {
  cout << endl << "Solve the system : A * x = b\n";
  cout << "en utisant  :\n";  
  cout <<"\t-Creation du  profil"<<endl;
  cout <<"\t-Algorithle Cuthill Mackee"<<endl;
  Cuthill cuthill = Cuthill("data.txt");
  int vertise = cuthill.getFirstNoeud();
  cuthill.cuthillMckceInverse(vertise);
  cuthill.optimizeProfil();
  cuthill.factorize();
  cuthill.resolutionInf();
  cuthill.resolutionSup();
  cuthill.solve();
  cuthill.displayResult();
  return 0;
}
