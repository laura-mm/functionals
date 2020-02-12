// changed z2 to match the paper
// also simplified parameters
// z1 <= z2
#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include <iomanip>
//#include <Eigen/Dense>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;
string carryon;
double a;
double mu;
bool info = false;

// if sig = 0 then z12 are infinite, and both param3 and param4 are infinite
double param3(double z1, double z2) // help/(sig*sqrt{q})
{
	if (a >= 1.0) return (z2 - z1)/(a + 1.0);
	else return (z2 - z1)/(2.0*a);
}
// if param3 = 0 then z1 = z2 and there is no middle section, everything is saturated, this happens when sig = inf
double param4(double z1, double z2) // (1+muM)/(sig*sqrt{q})
{
	if (a >= 1.0) return -z1;
	else return (((1.0 - a)*z2) - ((1.0 + a)*z1))/(2.0*a);
}
// if param4 = 0 then I suppose that means 1+muM = 0, or sig = inf
double M(double z1, double z2) // M
{
	//if (z1 == 1.0/0.0 && z2 == -1.0/0.0) return 1.0/(1.0 - mu); // help with mu
	//double p4 = param4(z1, z2);
	//if (p4 == 0.0) return -1.0/mu; // help with mu
	double m = (a + 1.0)*(1 - erf(z2/sqrt(2.0)))/2.0; // top saturation
	if (z1 != z2)
	{
		double p3 = param3(z1, z2);
		double p4 = param4(z1, z2);
		m += (p4*(erf(z2/sqrt(2.0)) - erf(z1/sqrt(2.0))))/(p3*2.0); // constant saturation
		m += (exp(-z1*z1/2.0) - exp(-z2*z2/2.0))/(p3*sqrt(2.0*M_PI)); // linear saturation
	}
	if (a < 1.0) m += (1.0 - a)*(1 + erf(z1/sqrt(2.0)))/2.0; // bottom saturation
	return m;
}
// M = 0 for a >= 1 and z2 = z1 = inf??
// what if mu = 1??
double q(double z1, double z2) // q
{
	//if (z1 == 1.0/0.0 && z2 == -1.0/0.0) return 1.0/pow(1.0 - mu, 2.0);
	double qu = pow(a + 1.0, 2.0)*(1.0 - erf(z2/sqrt(2.0)))/2.0;
	if (z1 != z2)
	{
		double p3 = param3(z1, z2); // sig*sqrt{q}
		double p4 = param4(z1, z2);
		qu += (pow(p4, 2.0) + 1.0)*(erf(z2/sqrt(2.0)) - erf(z1/sqrt(2.0)))/(2.0*pow(p3, 2.0)); // constant saturation
		qu += p4*(exp(-z1*z1/2.0) - exp(-z2*z2/2.0))/(pow(p3, 2.0)*sqrt(2.0*M_PI));
		qu -= (1.0 + a)*exp(-z2*z2/2.0)/(p3*sqrt(2.0*M_PI));
		if (a < 1.0) qu += (1.0 - a)*exp(-z1*z1/2.0)/(p3*sqrt(2.0*M_PI));
	}
	if (a < 1.0) qu += pow(1.0 - a, 2.0)*(1.0 - erf(z1/sqrt(2.0)))/2.0;
	return qu;
}
double sigma(double z1, double z2)
{
	if (q(z1, z2) == 0.0) return 0.0;
	else if (param4(z1, z2) == 0.0) return 1.0/0.0;
	else return (1.0 + mu*M(z1, z2))/(param4(z1, z2)*sqrt(q(z1, z2)));
}

double help(double z1, double z2)
{
	if (sigma(z1, z2) == 0.0) return 1.0;
	else return param3(z1, z2)*sigma(z1, z2)*sqrt(q(z1, z2));
}
double X(double z1, double z2)
{
	return (erf(z2/sqrt(2.0)) - erf(z1/sqrt(2.0)))/(2.0*help(z1, z2));
}
double gamma(double z1, double z2)
{
	return (1.0 - help(z1, z2))/(X(z1, z2)*pow(sigma(z1, z2), 2.0));
}
double phi(double z1, double z2)
{
	return (erf(z2/sqrt(2.0)) - erf(z1/sqrt(2.0)))/2.0;
}



//so algorithm to map z1 z2 to gamma and condition help = 1/(1+gamma)
//so in this algorithm we set mu, a, gamma
//and minimise gam - gamma and help - 1/(1+gamma) or 1/help - (1 + gamma)
// i chose to use help instead of chi and sigma, as help is required to find chi
//I think 2 newton raphsons at the same time maybe


// loop over z2, take many intervals
// for this z2 find the z1 that fufills stability criteria
// then find what gamma and sigma are
// and save these values
// order according to gamma
// i have simplified the condition to find phi = p3^2q


double gradient(double z1, double z2, double small) // with varying z1
{
	double grad = pow(param3(z1 + small, z2), 2.0)*q(z1 + small, z2);
	grad -= phi(z1 + small, z2);
	grad -= pow(param3(z1 - small, z2), 2.0)*q(z1 - small, z2);
	grad += phi(z1 - small, z2);
	return grad/(2.0*small);
}

double find_z1(double z1start, double z2)
{
	if (info == true)
	{
	ofstream zzz; zzz.open("zz.txt");
	ofstream ab; ab.open("abs.txt");

	cout << "now plotting function for findz1 for z2 = " << z2 << endl;
	
	for (int i = 0; i <= 200000; i++)
	{
		double zz = (0.0001*(double)i) - 10.0;
		zzz << zz;
		ab << (pow(param3(zz, z2), 2.0)*q(zz, z2)) - phi(zz, z2);
		if (i != 200000) {zzz << ", "; ab << ", ";}
	}
	zzz.close(); ab.close();
	cout << "finished plotting findz1" << endl;
	}
	
	double z1 = z1start; // this start could be from previous mu, not previous z2
	double small = pow(10.0, -6.0);
	double y = (pow(param3(z1, z2), 2.0)*q(z1, z2)) - phi(z1, z2);
	double grad;
	if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
	double ynew;
	double znew;
	while (abs(y) > 0.0)
	{
		grad = gradient(z1, z2, small);
		if (grad == 0.0) break; //{cout << "grad 0 " << g(z1, zz + small) << ", " << g(z1, zz - small) << endl;}
		if (info == true) {cout << "small " << small << " grad " << grad << endl;}
		znew = z1 - (y/grad);
		// if (z1 + zz < 0.0) cout << "badd" << endl;
		ynew = (pow(param3(znew, z2), 2.0)*q(znew, z2)) - phi(znew, z2);
		if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) break;
		y = ynew;
		z1 = znew;
		if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
		if (small > pow(10.0, -10.0)) small /= 10.0;
	}
	
	if (!(abs(sigma(z1, z2)) >= 0.0  && abs(gamma(z1, z2)) >= 0.0)) // found the wrong solution
	{
		small = pow(10.0, -6.0);
		grad = gradient(z1, z2, small);
		//z1 -= 2.0*y/grad;
		
		z1 *= -1.0;
		y = (pow(param3(z1, z2), 2.0)*q(z1, z2)) - phi(z1, z2);
		if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
		
		while (abs(y) > 0.0) // assume there are only 2 roots
		{
			grad = gradient(z1, z2, small);
			if (grad == 0.0) break; //{cout << "grad 0 " << g(z1, zz + small) << ", " << g(z1, zz - small) << endl;}
			if (info == true) {cout << "small " << small << " grad " << grad << endl;}
			znew = z1 - (y/grad);
			// if (z1 + zz < 0.0) cout << "badd" << endl;
			ynew = (pow(param3(znew, z2), 2.0)*q(znew, z2)) - phi(znew, z2);
			if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) break;
			y = ynew;
			z1 = znew;
			if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
			if (small > pow(10.0, -10.0)) small /= 10.0;
		}
	}
	
	return z1;	
}


bool compare(vector<double> i, vector<double> j) {return (i[0] < j[0]);} // only for 2 dimensional
//bool compare(ArrayXd i, ArrayXd j) {return (i(0) < j(0));} // only for 2 dimensional


void crit(int z2grid)
{
	vector<vector<double>> gamsig;
	//ArrayX2d gamsig;
	double start;
	for (int i = 0; i < z2grid; i ++)
	{
		//double z2 = exp((5.0*(double)i/(double)z2grid) - 3.0) - (mu);
		double z2 = (20.0*(double)i/(double)z2grid) - 10.0;
		//using e^x as is inverse of functio height and keeps change in area sort of constant
		// will need to see what kind of gamma this gives
		double z1start;
		if (i == 0) z1start = z2/a; // think about where this comes from
		else z1start = start;
		double z1 = find_z1(z1start, z2);
		start = z1;
		double g = gamma(z1, z2);
		double s = sigma(z1, z2);
		if (g >= -1.0 && g <= 1.0 && s >= 0)
		{
			cout << i << ", " << s << ", " << g << endl;
			vector<double> entry(2);
			entry[0] = g;
			entry[1] = s;
			gamsig.push_back(entry);
		}
	}
	
	sort(gamsig.begin(), gamsig.end(), compare);
	
	ofstream file;
	file.open("gamsig.txt");
	
	for (int j = 0; j < gamsig.size(); j ++) //cout << gamsig[j][0] << ", " << gamsig[j][1] << endl;
	{
		file << gamsig[j][0] << "," << gamsig[j][1];
		if (j != gamsig.size() - 1) file << ",";
	}
	file.close();
}



void m_crit(int z2grid, int mgrid, int ggrid)
{
	string filename = "m10crit_" + to_string((int)(10*a)) + ".txt";
	ofstream file; file.open(filename);	

	
	for (int m = 0; m <= mgrid; m++)
	{
		mu = (12.0*(double)m/(double)mgrid) - 2.0;
		
		vector<vector<double>> gamsig;
		double start;
		for (int iz = 0; iz < z2grid; iz ++)
		{
			//double z2 = exp((5.0*(double)iz/(double)z2grid) - 3.0) - (mu);
			double z2 = (10.0*(double)iz/(double)z2grid) - 0.0;
			//using e^x as is inverse of functio height and keeps change in area sort of constant
			// will need to see what kind of gamma this gives
			double z1start;
			if (iz == 0) z1start = z2/a; // think about where this comes from
			else z1start = start;
			double z1 = find_z1(z1start, z2);
			start = z1;
			double g = gamma(z1, z2);
			double s = sigma(z1, z2);
			if (g >= -1.0 && g <= 1.0 && s >= 0)
			{
				//cout << i << ", " << s << ", " << g << endl;
				vector<double> entry(2);
				entry[0] = g;
				entry[1] = s;
				gamsig.push_back(entry);
			}
		}
	
		sort(gamsig.begin(), gamsig.end(), compare);

		for (int ig = 1; ig <= ggrid; ig ++)
		{
			cout << m << ", " << ig << endl;
			
			double gam = (2.0*(double)ig/(double)ggrid) - 1.0;
			ArrayXd absol(gamsig.size());
			for (int j = 0; j < gamsig.size(); j++) absol(j) = abs(gamsig[j][0] - gam);
			//vector<double> absol(gamsig.size());
			//for (int j = 0; j < gamsig.size(); j++) absol[j] = abs(gamsig[j][0] - gam);
			int location;
			double min = absol.minCoeff(&location);
			if (min > pow(10.0, -2.0)) {cout << "mu = " << mu << " gam = " << gam << " min = " << min << endl; exit(1);}
			else
			{
				file << gamsig[location][1];
				if (ig != ggrid) file << ",";
			}
		}
		if (m != mgrid) file << ",";
	}
	file.close();
}












int main()
{
	double zgrid = 100000; // number of z2's to try
	double mgrid = 200;
	double ggrid = 200;
	a = 0.5;
	mu = -0.0;
	
	crit(zgrid);
	
	//m_crit(z2grid, mgrid, ggrid);
	/*
	double z2 = 1.0;
	double z1start = 1.0;
	double z1 = find_z1(z1start, z2);
	cout << "z1 = " << z1 << " z2 = " << z2 << " gamma = " << gamma(z1, z2) << " sigma = " << sigma(z1, z2) << endl;
	*/


	
	return 0;
}






