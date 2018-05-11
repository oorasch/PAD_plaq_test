#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <string>
#include <math.h>
#include <limits>
#include <gsl/gsl_sf_bessel.h>

using namespace std;

const int V = 16;
const int max_bessel = 50;

void bessel_init(const double& eta, const double& eta_bar, vector<double>& vec);
double meas_plaq(const double& eta, const double& eta_bar, const vector<int>& plaq_occ, const vector<double>& I_bessel);

int main(){

int range = 0;
int q = 3;
vector<int> plaq_occ;
vector<double> I_bessel;

for(int i = 0; i < V; i++) plaq_occ.push_back(q);

 stringstream s1;
 s1 << q;
 string s_q = s1.str();

ofstream outfile1, outfile2;
outfile1.open("data/plaq_test_"+s_q+".dat");

for(int i = 1; i <= 100; i++)
{
	double eta = 0.5*i*0.1, eta_bar = 0.5*i*0.1;
	
	bessel_init(eta, eta_bar, I_bessel);
	
	//for(int i= -1 ; i < 1; i++)
	{
		//vector<int> plaq_occ (V, i);
		outfile1 << 2.0*eta << " " << meas_plaq(eta, eta_bar, plaq_occ, I_bessel) << endl;
	
	//cout << "#############################################" << endl;
	
	}
	
}
outfile1.close();

return 0;

}

void bessel_init(const double& eta, const double& eta_bar, vector<double>& vec)
{
	vec.clear();
	double aux = 2.0*sqrt(eta*eta_bar);
	//cout << aux << endl;
	for(int i = 0; i < max_bessel; i++)
	{
		vec.push_back(gsl_sf_bessel_In(i, aux));
		//cout << vec.at(i) << " ";
	}
	//cout << endl;
}

double meas_plaq(const double& eta, const double& eta_bar, const vector<int>& plaq_occ, const vector<double>& I_bessel)
{
	double tmp = 0.0;
	double aux1 = 0.25*(eta + eta_bar)/sqrt(eta*eta_bar),
				 aux2 = (0.25/eta - 0.25/eta_bar);

	for(auto it: plaq_occ)
	{
		if(it == 0)
		{
			tmp += 2.0*aux1*I_bessel.at(1)/I_bessel.at(0);
		}
		else
		{
			tmp += aux1*(I_bessel.at(abs(it)+1) + I_bessel.at(abs(it)-1))/I_bessel.at(abs(it));
			tmp += it*aux2;
		} 				
		
	}
	
	return tmp/V;	
}
