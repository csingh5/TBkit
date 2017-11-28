///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** main.cpp ***                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include "Eigmesh.hpp"
#include "Mesh.hpp"
#include "W90.hpp"

using namespace std;

inline double fermi(double epsilon) {
	return 1.0 / (exp(1000*epsilon)   + 1.0);
}

int main(int argc, char *argv[])
{
	
	Eigmesh eigmesh("data/SmB6.klist");
	W90 w90("data/wan_hr.dat");
	
	for (int i = 0; i < eigmesh.wedge(); ++i) {
		saes_t s(w90.fft(eigmesh.point(i)));
		eigmesh.push_back(s);	
	}

	Eigen::ArrayXd d = Eigen::ArrayXd::Zero(18);

	for (int i = 0; i < eigmesh.wedge(); ++i) {
		for (int band = 0; band < 18; ++band) {
			d(band) += eigmesh.weight(i) * fermi(eigmesh.eigenvalue(i, band));
		}
	}

	d /= eigmesh.full();
	cout << d.sum() << endl;

	d *= 0;

	for (int i = 0; i < eigmesh.wedge(); ++i) {
		for (int band = 0; band < 18; ++band) {
			for (int vec = 0; vec < 18; ++vec) {
				d(band) += 
					eigmesh.weight(i) * 
					fermi(eigmesh.eigenvalue(i, vec)) *
					eigmesh.abs2eigenvector(i, band, vec);
			}
		}
	}
	
	d /= eigmesh.full();

	cout << d.head(9) + d.tail(9) << endl;
	cout << d.sum() << endl;



	return 0;
}
