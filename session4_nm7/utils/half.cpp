/*
     -------------------------------------------------------------
        Program to calculate half-life for two-compartment model
     -------------------------------------------------------------
*/


#include <iostream>
#include <cmath>


using namespace std;

int main() {
        double cl,vc,q,vp; 

cout << "To obtain alpha and beta half-lifes enter the values" << endl;
cin >> cl; 
cin >> vc; 
cin >> q; 
cin >> vp;

double kel, k12, k21, val1,val2, alpha,beta; 

        kel = cl/vc;
        k12 = q/vc;
        k21 = q/vp;
        
        val1 = kel + k12 + k21;
        val2 = kel * k21; 

        alpha = (val1 + sqrt(pow(val1,2) - 4 * val2))/2;
        beta  = (val1 - sqrt(pow(val1,2) - 4 * val2))/2;
        double hla = 0.693/alpha;
        double hlb = 0.693/beta; 

        cout << "Alpha half life :" << hla << endl;
        cout << "Beta half life :" << hlb  << endl;

return 0;
}

