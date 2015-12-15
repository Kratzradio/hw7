#include<cmath>
#include<iostream>
#include<fstream>

using namespace std;

void ssc(const double* f, const double* rk5, double& max);

void RK45(double* f, double* rk5, const double mu, double dt, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7);

void RK(const double* f, const double mu, const double dt, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7);

void function(double* const f1, const double mu, const double x, const double y, const double x1, const double y1);

int main(){
ofstream out("solution");
const int dim = 4;
double dt = 1e-4;
double T = 17.065216560157;
double t = 0.0;
double mu = 0.012277471;
double f[dim] = {0.994, 0.0, 0.0, -2.00158510637908};
double rk5[dim] = {0.994, 0.0, 0.0, -2.00158510637908};
double k1[dim], k2[dim], k3[dim], k4[dim], k5[dim], k6[dim], k7[dim];
double max;
double tol = 1e-5;
out << t << "\t" << f[0] << "\t" << f[1] << "\t" << rk5[0] << "\t" << rk5[1] << "\t" << max << "\t" << dt << endl;
	while(t<=T){
	RK45(f, rk5, mu, dt, k1, k2, k3, k4, k5, k6, k7);
	ssc(f, rk5, max);
	dt*=pow(tol/max, 0.2);	
	t+=dt;
out << t << "\t" << f[0] << "\t" << f[1] << "\t" << rk5[0] << "\t" << rk5[1] << "\t" << max << "\t" << dt << endl;
	}

out.close();
return 0;
}

void ssc(const double* f, const double* rk5, double& max){
double e = 0.0;
max = 0.0;
for(int i=0; i<4; i++){
	e = abs(rk5[i]-f[i]);
	if(e>max) max = e;
	}

}

void RK45(double* f, double* rk5, const double mu, double dt, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7){

const double b41=5179.0/57600.0, b42=0.0, b43=7571.0/16695.0, b44=393.0/640.0, b45=-92097.0/339200.0, b46=187.0/2100.0, b47=1.0/40.0;
const double b51=35.0/384.0, b52=0.0, b53=500.0/1113.0, b54=125.0/192.0, b55=-2187.0/6784.0, b56=11.0/84.0;
RK(f, mu, dt, k1, k2, k3, k4, k5, k6, k7);
for(int i=0; i<4; i++){
	rk5[i] = f[i] + dt*(b51*k1[i]+b52*k2[i]+b53*k3[i]+b54*k4[i]+b55*k5[i]+b56*k6[i]);
	f[i] += dt*(b41*k1[i]+b42*k2[i]+b43*k3[i]+b44*k4[i]+b45*k5[i]+b46*k6[i]+b47*k7[i]);
	}
}

void RK(const double* f, const double mu, const double dt, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7){

	double x=f[0], y=f[1], x1=f[2], y1=f[3];

	const double a21=1.0/5.0;
	const double a31=3.0/40.0, 	 a32=9.0/40.0;
	const double a41=44.0/45.0, 	 a42=-56.0/15.0,      a43=32.0/9.0;
	const double a51=19372.0/6561.0, a52=-25360.0/2187.0, a53=64448.0/6561.0, a54=-212.0/729.0;
	const double a61=9017.0/3168.0,  a62=-355.0/33.0,     a63=46732.0/5247.0, a64=49.0/176.0,  a65=-5103.0/18656.0;
	const double a71=35.0/384.0, 	 a72=0.0,	      a73=500.0/1113.0,   a74=125.0/192.0, a75=-2187.0/6784.0, a76=11.0/84.0;

	function(k1, mu, x, y, x1, y1);
	function(k2, mu, x+dt*(a21*k1[0]), y+dt*(a21*k1[1]), x1+dt*(a21*k1[2]), y1+dt*(a21*k1[3]));
	function(k3, mu, x+dt*(a31*k1[0]+a32*k2[0]), y+dt*(a31*k1[1]+a32*k2[1]), x1+dt*(a31*k1[2]+a32*k2[2]), y1+dt*(a31*k1[3]+a32*k2[3]));
	function(k4, mu, x+dt*(a41*k1[0]+a42*k2[0]+a43*k3[0]), y+dt*(a41*k1[1]+a42*k2[1]+a43*k3[1]), x1+dt*(a41*k1[2]+a42*k2[2]+a43*k3[2]), y1+dt*(a41*k1[3]+a42*k2[3]+a43*k3[3]));
	function(k5, mu, x+dt*(a51*k1[0]+a52*k2[0]+a53*k3[0]+a54*k4[0]), y+dt*(a51*k1[1]+a52*k2[1]+a53*k3[1]+a54*k4[1]), x1+dt*(a51*k1[2]+a52*k2[2]+a53*k3[2]+a54*k4[2]), y1+dt*(a51*k1[3]+a52*k2[3]+a53*k3[3]+a54*k4[3]));
	function(k6, mu, x+dt*(a61*k1[0]+a62*k2[0]+a63*k3[0]+a64*k4[0]+a65*k5[0]), y+dt*(a61*k1[1]+a62*k2[1]+a63*k3[1]+a64*k4[1]+a65*k5[1]), x1+dt*(a61*k1[2]+a62*k2[2]+a63*k3[2]+a64*k4[2]+a65*k5[2]), y1+dt*(a61*k1[3]+a62*k2[3]+a63*k3[3]+a64*k4[3]+a65*k5[3]));
	function(k7, mu, x+dt*(a71*k1[0]+a72*k2[0]+a73*k3[0]+a74*k4[0]+a75*k5[0]+a76*k6[0]), y+dt*(a71*k1[1]+a72*k2[1]+a73*k3[1]+a74*k4[1]+a75*k5[1]+a76*k6[1]), x1+dt*(a71*k1[2]+a72*k2[2]+a73*k3[2]+a74*k4[2]+a75*k5[2]+a76*k6[2]), y1+dt*(a71*k1[3]+a72*k2[3]+a73*k3[3]+a74*k4[3]+a75*k5[3]+a76*k6[3]));

}

void function(double* const f1, const double mu, const double x, const double y, const double x1, const double y1){
	double r = sqrt((x+mu)*(x+mu)+y*y);
	double s = sqrt((x-1.0+mu)*(x-1.0+mu)+y*y);
	f1[0] = x1;
	f1[1] = y1;
	f1[2] = x + 2.0*y1 - (1.0-mu)*(x+mu)/(r*r*r) - mu*(x-1.0+mu)/(s*s*s);
	f1[3] = y - 2.0*x1 - (1.0-mu)*y/(r*r*r) - mu*y/(s*s*s);


}
