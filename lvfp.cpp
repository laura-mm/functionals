// fixed point solver using trick method for lv
// this works for mu = 0 only!
// z1 is the same for both cases but z2 different
// in this program z1 >= z2, so opposite to the paper

#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include <eigen3/Eigen/Dense>
#include <iomanip>
using namespace Eigen;
using namespace std;

IOFormat p(20, DontAlignCols, ",", ",", "", "", "", "");
IOFormat c(6, DontAlignCols, ",", ",", "", "", "", "");

double a;
double mu;

bool info;



double param1(double z1, double z2) // (1+muM)/2(1-gXsig^2)
{
	if (a >= 1.0) {return (z2*(a + 1.0))/(2.0*(z1 + z2));}
	else {return (z1 + z2 + a*(z2 - z1))/(2.0*(z1 + z2));}
}
double param2(double z1, double z2) // (sig*sqrt{q})/sqrt{2pi}(1-gXsig^2)
{
	if (a >= 1.0) {return (a + 1.0)/(sqrt(2.0*M_PI)*(z1 + z2));}
	else {return (sqrt(2.0)*a)/(sqrt(M_PI)*(z1 + z2));}
}
double M(double z1, double z2) // M
{
	if (z1 == 1.0/0.0 || z2 == 1.0/0.0) return 1.0/(1.0 - mu);
	double m = (a + 1.0)*(1 - erf(z1/sqrt(2.0)))/2.0;
	m += param1(z1, z2)*(erf(z1/sqrt(2.0)) + erf(z2/sqrt(2.0)));
	m += param2(z1, z2)*(exp(-z2*z2/2.0) - exp(-z1*z1/2.0));
	if (a >= 1.0) {return m;}
	else {return m + ((1.0 - a)*(1 - erf(z2/sqrt(2.0)))/2.0);}
}
double help(double z1, double z2) // (1-gXsig^2)
{
	if (z1 + z2 == 0.0)
	{
		double p; // p2/p1
		if (a >= 1.0) p = sqrt(2.0)/(sqrt(M_PI)*z2);
		else p = sqrt(8.0)/(sqrt(M_PI)*(z2 - z1));
		double h = (mu*(erf(z1/sqrt(2.0)) + erf(z2/sqrt(2.0))))/2.0;
		h += (mu*p*(exp(-z2*z2/2.0) - exp(-z1*z1/2.0)))/2.0;
		return h;
	}
	else return (1.0 + (mu*M(z1, z2)))/(2.0*param1(z1, z2));
}
double var(double z1, double z2) // sig*sqrt{q}
{
	return param2(z1, z2)*sqrt(2.0*M_PI)*help(z1, z2);
}
double X(double z1, double z2) // X
{
	return (erf(z1/sqrt(2.0)) + erf(z2/sqrt(2.0)))/(2.0*help(z1, z2));
}
double q(double z1, double z2) // q
{
	if (z1 == 1.0/0.0 || z2 == 1.0/0.0) return 1.0/pow(1.0 - mu, 2.0);
	double qu = (pow(a + 1.0, 2.0)*(1.0 - erf(z1/sqrt(2.0)))/2.0);
	qu += (((2.0*pow(param1(z1, z2), 2.0)) + (M_PI*pow(param2(z1, z2), 2.0)))*(erf(z1/sqrt(2.0)) + erf(z2/sqrt(2.0))));
	qu += 2.0*param1(z1, z2)*param2(z1, z2)*(exp(-z2*z2/2.0) - exp(-z1*z1/2.0));
	qu -= (1.0 + a)*param2(z1, z2)*exp(-z1*z1/2.0);
	if (a >= 1.0) {return qu;}
	else
	{
		qu += (1.0 - a)*param2(z1, z2)*exp(-z2*z2/2.0);
		qu += pow(1.0 - a, 2.0)*(1.0 - erf(z2/sqrt(2.0)))/2.0;
		return qu;
	}
}
double sig(double z1, double z2) // sigma
{
	if (z1 == 1.0/0.0 || z2 == 1.0/0.0) return 0.0;
	else if (z1 + z2 == 0.0)
	{
		double p; // p1/p2
		if (a >= 1.0) p = (sqrt(M_PI)*z2)/sqrt(2.0);
		else p = (sqrt(M_PI)*(z2 - z1))/sqrt(8.0);
		double q = (pow(p, 2.0) + (M_PI/2.0))*(erf(z1/sqrt(2.0)) + erf(z2/sqrt(2.0)));
		q += p*(exp(-z2*z2/2.0) - exp(-z1*z1/2.0));
		q /= M_PI*pow(help(z1, z2), 2.0);
		return pow(q, -0.5);
	}
	else return var(z1, z2)/sqrt(q(z1, z2));
}
double g(double z1, double z2) // gamma
{
	return (1.0 - help(z1, z2))/(X(z1, z2)*pow(sig(z1, z2), 2.0));
}
double top(double z1) {return (1.0 - erf(z1/sqrt(2.0)))/2.0;}
double mid(double z1, double z2) {return (erf(z1/sqrt(2.0)) + erf(z2/sqrt(2.0)))/2.0;}
double bot(double z1, double z2) {return (1.0 - erf(z2/sqrt(2.0)))/2.0;}



double find_z2(double z1, double gam) // without a start
{
	/*ofstream zzz; zzz.open("zz.txt");
	ofstream ab; ab.open("abs.txt");

	for (int i = 0; i <= 20000; i++)
	{
		double zz = (0.001*(double)i) - 10.0;
		zzz << zz;
		//cout << endl << zz << endl;
		ab << g(z1, zz, a, mu) - gam;
		// cout << "zz " << zz << " ymax " << ymax << endl;
		if (i != 20000) {zzz << ", "; ab << ", ";}
	}*/

	double k = (((1.0 - a)*(1.0 - mu)) + (2.0*mu))/((1.0 - mu)*(1.0 + a));
	double zz = (1.0 + k)*z1/(1.0 - k);
	//if (z1 + zz < 0.0) cout << "badd" << endl;
	double small = 0.1;
	double y = g(z1, zz) - gam;
	double grad;
	//cout << "zz " << zz << ", y " << y << endl;;
	double ynew;
	while (abs(y) > 0.0)
	{
		grad = (g(z1, zz + small) - g(z1, zz - small))/(2.0*small);
		//cout << "grad " << grad << endl;
		zz -= (y/grad);
		// if (z1 + zz < 0.0) cout << "badd" << endl;
		ynew = g(z1, zz) - gam;
		if (abs(ynew) >= abs(y) && abs(y) < pow(10.0, -10.0)) break;
		y = ynew;
		//cout << "zz " << zz << ", y " << y << endl;
		small /= 1.5;
	}
	return zz;
		
}
double find_z2(double z1, double gam, double start) // with a start
{

	if (info == true) {cout << "function find z2" << endl;}

	if (info == true)
	{
	ofstream zzz; zzz.open("zz.txt");
	ofstream ab; ab.open("abs.txt");

	for (int i = 0; i <= 20000000; i++)
	{
		double zz = (0.000001*(double)i) - 10.0;
		zzz << zz;
		//cout << endl << zz << endl;
		ab << g(z1, zz) - gam;
		if (zz == 0.0)
		{
			cout << "here gamma = " << g(z1, zz) << endl;
		}
		// cout << "zz " << zz << " ymax " << ymax << endl;
		if (i != 20000000) {zzz << ", "; ab << ", ";}
	}
	}

	//double k = (((1.0 - a)*(1.0 - mu)) + (2.0*mu))/((1.0 - mu)*(1.0 + a));
	double zz = start; // (1.0 + k)*z1/(1.0 - k);
	//if (z1 + zz < 0.0) cout << "badd" << endl;
	double small = pow(10.0, -7.0);
	double y = g(z1, zz) - gam;
	double grad;
	if (info == true) {cout << "zz " << zz << ", y " << y << endl;}
	double ynew;
	while (abs(y) > 0.0)
	{
		grad = (g(z1, zz + small) - g(z1, zz - small))/(2.0*small);
		if (grad == 0.0 && abs(y) < pow(10.0, -5.0)) break; //{cout << "grad 0 " << g(z1, zz + small) << ", " << g(z1, zz - small) << endl;}
		if (info == true) {cout << "small " << small << " grad " << grad << endl;}
		zz -= (y/grad);
		// if (z1 + zz < 0.0) cout << "badd" << endl;
		ynew = g(z1, zz) - gam;
		if (abs(ynew) >= abs(y) && abs(y) < pow(10.0, -10.0)) break;
		y = ynew;
		if (info == true) {cout << "zz " << zz << ", y " << y << endl;}
		if (small > pow(10.0, -10.0)) small /= 10.0;
	}
	return zz;
	if (info == true) {cout << "function find z2 end" << endl;}
		
}
vector<double> yy(double z1, double gam) // [0] yy
{
	double z2 = find_z2(z1, gam);
	//cout << z2 << endl;
	//cout << X(z1, z2, a, mu) << endl;
	//cout << sig(z1, z2, a, mu) << endl;
	vector<double> vec(1);
	vec[0] = (pow(X(z1, z2), -1.0)*pow(sig(z1, z2), -2.0)) - (1.0 + gam);
	return vec;
}
vector<double> yy(double z1, double gam, double start_z2) // [0] yy, [1] z2
{
	double z2 = find_z2(z1, gam, start_z2);
	if (info == true) {cout << z1 << ", " << start_z2 << ", " << z2 << endl;}
	/*if (param2(z1, z2, a) < 0.0)
	{
		cout << "try the other one" << endl;
		double newstart;
		if (a >= 1.0) newstart = -z1*z2/(z1 + (2.0*z2));
		else newstart = z1*(((a-1.0)*z1) - z2)/(z1 + ((1.0 + a)*z2));
		z2 = find_z2(z1, gam, a, mu, newstart);
	}*/

	//cout << z2 << endl;
	//cout << X(z1, z2, a, mu) << endl;
	//cout << sig(z1, z2, a, mu) << endl;
	vector<double> vec(2);
	vec[0] = (pow(X(z1, z2), -1.0)*pow(sig(z1, z2), -2.0)) - (1.0 + gam);
	vec[1] = z2;
	return vec;
}
vector<double> crit(double gam) // [0] sig, [1] z1 
{
	/*ofstream zzz; zzz.open("zz.txt");
	ofstream other; other.open("other.txt");
	ofstream ab; ab.open("abs.txt");

	for (int i = 0; i <= 20000; i++)
	{
		double zz = (0.001*(double)i) - 10.0;
		zzz << zz;
		other << find_z2(zz, gam, a, mu);
		//cout << endl << zz << endl;
		ab << yy(zz, gam, a, mu)[0];
		// cout << "zz " << zz << " ymax " << ymax << endl;
		if (i != 20000) {zzz << ", "; ab << ", "; other << ", ";}
	}*/
	double zz = 1.0;
	double small = 0.1;
	vector<double> vec = yy(zz, gam);
	double y = vec[0];
	double grad;
	//cout << "zz " << zz << ", y " << y << endl;;
	double ynew;
	double znew;
	while (abs(y) > 0.0)
	{
		grad = (yy(zz + small, gam)[0] - yy(zz - small, gam)[0])/(2.0*small);
		//cout << "grad " << grad << endl;
		znew = zz - (y/grad);
		//zz -= (y/grad);
		vec = yy(znew, gam);
		ynew = vec[0];
		if (abs(ynew) >= abs(y) || !(abs(ynew) >= 0.0)) break;
		zz = znew;
		y = ynew;
		//cout << "zz " << zz << ", y " << y << endl;
		small /= 1.5;
	}
	vector<double> vecc(3);
	vecc[0] = sig(zz, find_z2(zz, gam));
	vecc[1] = zz;
	vecc[2] = find_z2(zz, gam);
	return vecc;
}
vector<double> crit(double gam, double start_z1, double start_z2) // [0] sig, [1] z1, [2] z2
{
	/*ofstream zzz; zzz.open("zz.txt");
	ofstream other; other.open("other.txt");
	ofstream ab; ab.open("abs.txt");

	for (int i = 0; i <= 20000; i++)
	{
		double zz = (0.001*(double)i) - 10.0;
		zzz << zz;
		other << find_z2(zz, gam, a, mu, start_z2);
		//cout << endl << zz << endl;
		ab << yy(zz, gam, a, mu, start_z2)[0];
		// cout << "zz " << zz << " ymax " << ymax << endl;
		if (i != 20000) {zzz << ", "; ab << ", "; other << ", ";}
	}*/
	double zz = start_z1;
	double small = 0.0001;
	vector<double> vec = yy(zz, gam, start_z2);
	double y = vec[0];
	double z2 = vec[1];
	double grad;
	if (info == true) {cout << "zz " << zz << ", y " << y << endl;}
	double ynew;
	double znew;
	double z2new;
	while (abs(y) > 0.0)
	{
		grad = (yy(zz + small, gam, z2)[0] - yy(zz - small, gam, z2)[0])/(2.0*small);
		if (info == true) {cout << "grad " << grad << endl;}
		znew = zz - (y/grad);
		//zz -= (y/grad);
		vec = yy(znew, gam, z2);
		ynew = vec[0];
		z2new = vec[1];
		if (abs(ynew) >= abs(y) || !(abs(ynew) >= 0.0)) break;
		zz = znew;
		y = ynew;
		z2 = z2new;
		if (info == true) {cout << "zz " << zz << ", y " << y << endl;}
		small /= 1.5;
	}
	vector<double> vecc(3);
	vecc[0] = sig(zz, z2);
	vecc[1] = zz;
	vecc[2] = z2;
	return vecc;
}



bool compare(ArrayXd i, ArrayXd j) {return (i(0) < j(0));} // cheack all this

void fiveplot(int grid) // now does z1 and z2
{
	string critstr = "crit_" + to_string((int)(10*a)) + "_" + to_string((int)(10*mu)) + ".txt";
	ofstream critfile; critfile.open(critstr);

	for (int gi = 0; gi <= 4; gi++)
	{
		string filename = "meashist_" + to_string((int)(10*a)) + "_" + to_string((int)(10*mu)) + "_" + to_string(gi) + ".txt"; // changed to hist for sigma != 0
		ofstream file; file.open(filename);


		double gamma = (0.5*(double)gi) - 1.0;
		critfile << crit(gamma)[0];

		vector<ArrayXd> sigmat;
		double start = -40.0;
		double end = 40.0;
		double range = end - start;

		for (int i = 0; i <= grid; i++)
		{
			double z1 = (a - mu)/(((double)i*range/(double)grid) + start);
			double z2 = find_z2(z1, gamma);
			if (abs(g(z1, z2) - gamma) > pow(10.0, -10.0)) continue; //{cout << "bad solution1 " << z1 << endl; continue;}
			if (sig(z1, z2) < 0.0) continue; //{cout << "badness" << endl; continue;}
			if (abs(mid(z1, z2) - help(z1, z2)*X(z1, z2)) > pow(10.0, -10.0)) continue; //{cout << "bad solution3 " << z1 << endl; continue;}
			if (!(abs(sig(z1, z2)) >= 0.0)) continue;
			if (sig(z1, z2) == 0.0) continue; // this is for histograms only

			ArrayXd entry = ArrayXd::Zero(9); // sigma, z1, z2, top, middle, bottom, M, q, help
			entry(0) = sig(z1, z2);
			entry(1) = z1;
			entry(2) = z2;
			entry(3) = top(z1);
			entry(4) = mid(z1, z2);
			entry(5) = bot(z1, z2);
			entry(6) = M(z1, z2);
			entry(7) = q(z1, z2);
			entry(8) = help(z1, z2);
			sigmat.push_back(entry);
			//cout << "i " << i << ", z1 " << z1 << ", z2 " << z2 << ", sig " << sig(z1, z2, a, mu) << ", 1/q " << 1.0/q(z1, z2, a, mu) << endl;
		}
		sort(sigmat.begin(), sigmat.end(), compare);

		cout << gi << " " << sigmat.size() << " sorted!" << endl;
	
		for (int i = 0; i < sigmat.size(); i++)
		{
			file << sigmat[i].format(c);
			//cout << sigmat[i].format(c) << endl;
			if (i != sigmat.size()-1) file << ",";
		}
		file.close();
		if (gi != 4) critfile << ",";
	}
	critfile.close();
}
void a_crit(int gridg, int grida)
{
	string filename = "acrit_" + to_string((int)(10*mu)) + ".txt";
	ofstream file; file.open(filename);
	for (int gi = 3; gi <= gridg; gi++)
	{
		double gamma =  (2.0*gi/(double)gridg) - 1.0;
		for (int ai = 3; ai <= grida; ai++)
		{
			a = 2.0*ai/(double)grida;
			vector<double> v = crit(gamma);
			if (v[0] <= 0.0) cout << "badness g = " << gi << ", a = " << ai << endl;
			file << v[0];
			if (ai != grida) file << ",";
		}
		if (gi != gridg) file << ",";
	}
	file.close();
}



void phase_line(int gridg)
{
	string filename = "phline_" + to_string((int)(10*a)) + "_" + to_string((int)(10*mu)) + ".txt";
	ofstream file; file.open(filename);
	for (int ig = 1; ig <= gridg; ig++)
	{
		double gamma = (2.0*(double)ig/(double)gridg) - 1.0;
		file << crit(gamma)[0];
		if (ig != gridg) file << ",";
	}
	file.close();
}



		

int main()
{



	int grid = 20; //for phase suprhace
	int grids = 400000; // for graphs with varying sigma
	int grida = 200; // for acrit, a grid
	int gridg = 200; // for acrit gamma grid
	//double gam = 0.0;
	a = 2.0;
	//mu = 0.0;

	//fiveplot(grids);


	//a_crit(gridg, gridg);

	m_crit(200, 40);

	//phase_line(1000);

	//mu = 0.435;
	//vector<double> v = crit(0.0, z1bad, z2bad);
	//cout << endl<< v[0] << ", " << v[1] << ", " << v[2] << endl;
	//cout << z1bad << ", " << z2bad << endl;



	return 0;
}


	/*double z2 = find_z2(z1, gam, a, mu, -0.205);
	//double z2 = -1.0;
	cout << "z1 " << z1 << " z2 " << z2 << endl;
	cout << "p1 = " << param1(z1, z2, a) << endl;
	cout << "p2 = " << param2(z1, z2, a) << endl;
	cout << "M = " << M(z1, z2, a, mu) << endl;
	cout << "help = " << help(z1, z2, a, mu) << endl;
	cout << "var = " << var(z1, z2, a, mu) << endl;
	cout << "X = " << X(z1, z2, a, mu) << endl;
	cout << "q = " << q(z1, z2, a, mu) << endl; // this one is bad
	cout << "sig = " << sig(z1, z2, a, mu) << endl;
	cout << "gam = " << g(z1, z2, a, mu) << endl;
	cout << "thing = " << 1.0 + (mu*M(z1, z2, a, mu)) << endl;


	cout << log10(sig(z1, z2, a, mu)) << endl;
*/
































		




