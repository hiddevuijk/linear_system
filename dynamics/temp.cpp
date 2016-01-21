// nr3 headers
#include "nr_headers/nr3.h"
#include "nr_headers/ran.h"
#include "nr_headers/fourier.h"
#include "nr_headers/correl.h"
#include "nr_headers/stepper.h"
#include "nr_headers/stepperbs.h"
#include "nr_headers/odeint.h"
#include "nr_headers/stepperdopr853.h"
#include "nr_headers/stepperdopr5.h"
#include "nr_headers/stepperstoerm.h"

// other headers
#include "headers/derivatives.h"
#include "headers/write_matrix.h"

// std headers
#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
using namespace::std;

int main(int argc, char *argv[])
{
	int N =50;

	int dt = 1000;
	if(argc == 2) {
		try {
			dt = stoi(argv[1]);
		}
		catch(exception& e) {
			cerr << "Unable to cast " << argv[1] << " to int" << endl;
			cerr << "ERROR: " << e.what() << endl;
			return 1;
		}
	}
		
	double ti = 10000;
	double tf = ti+dt;

	double th = 2.5;

	VecDoub xw(N+N*N);
	Ran r(12345677);
	for(int i=0;i<(N+N*N);++i)  xw[i] = r.doub();
	double alpha = 1e-3;	
	NW nw(alpha,N);

	Output out_0_ti;	
	Output out(dt);

	double atol = 1.e-5;
	double rtol = atol;
	double h1 = 0.01;
	double hmin = 0.0;

	// integrate from 0 to ti
	Odeint<StepperDopr5<NW> > ode_0_ti(xw,0,ti,atol,rtol,h1,hmin,out_0_ti,nw);
	ode_0_ti.integrate();

	// integrate from ti to tf
	Odeint<StepperDopr5<NW> > ode(xw,ti,tf,atol,rtol,h1,hmin,out,nw);
	ode.integrate();

	vector<int> sau(dt,0);
	vector<int> avs(dt,0);
	int j=0;
	for(int t=0;t<dt;++t){
		for(int i=0;i<N;++i) {
			if( abs(out.ysave[i][t]) > th){
				++sau[t];
			}
		}
		if(sau[t] != 0){
			avs[j] += sau[t];
		} else if ( avs[j] != 0) {
			++j;
		}
	}

	write_matrix(sau,dt,"sau.csv");
	write_matrix(avs,j,"avs.csv");
//	write_matrix(out.xsave,dt,"t.csv");
	write_matrix(out.ysave,N+N*N,dt,"xw.csv");

	return 0;
}
