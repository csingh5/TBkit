///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** w90.hpp ***                                  //
//                                                                           //
// class to manage methods relating to wannier90 hamiltonian                 //
// templatized over dimension of the mesh                                    //
//                                                                           //
// created November 27, 2017                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef W90_hpp
#define W90_hpp

#include <cassert>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

class W90
{
public:
	W90(std::string hr_file) {
		
		std::ifstream file(hr_file);
		assert(file.is_open());

		std::string line;
		double data;
		std::vector<double> v;

		while (std::getline(file, line)) {
			std::stringstream input(line);
			for (int i = 0; i < 7; ++i) {
				input >> data;
				v.push_back(data);
			}
		}

		hr = Eigen::Map<Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, 
		   Eigen::RowMajor> >(&v[0], v.size()/7, 7);
	}

	Eigen::MatrixXcd fft(const Eigen::Vector3d& k) {

		int dimension = hr.col(3).maxCoeff();
		Eigen::MatrixXcd h = Eigen::MatrixXcd::Zero(dimension, dimension);

		for (int nn = 0; nn < hr.rows(); ++nn) {
			h((int)hr(nn, 3)-1, (int)hr(nn, 4)-1) += tij(nn, k);
		}

		return h;
	}

private:
	Eigen::ArrayXXd hr;
	
	inline std::complex<double> tij(int nn, const Eigen::Vector3d& k) {
		
		std::complex<double> i(0, 1);
		std::complex<double> hopping(hr(nn, 5), hr(nn, 6));
		std::complex<double> phase;

		phase = exp(i * (k(0)*hr(nn, 0) + k(1)*hr(nn,1) + k(2)*hr(nn, 2)));

		return hopping * phase;	
	}
};


#endif /* W90_hpp */
