///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** Eigmesh.hpp ***                              //
//                                                                           //
// class to hold the eigenvalues and eigenvectors on a mesh                  //
//                                                                           //
// created November 27, 2017                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Eigmesh_hpp
#define Eigmesh_hpp

#include "Mesh.hpp"

using saes_t = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>;
using eV = double;

class Eigmesh : public Mesh
{
public:
	Eigmesh(std::string klist) : Mesh(klist) {
		m_eig.reserve(size());
	}
	
	eV eigenvalue(int mesh_point_index, int band_index) {
		assert(m_eig.size() != 0 && "Run a diagonalization first!");
		return m_eig[mesh_point_index].eigenvalues()(band_index);	
	}	

	std::complex<double> eigenvector(int mesh_point_index, int row, int col) {
		assert(m_eig.size() != 0 && "Run a diagonalization first!");
		return m_eig[mesh_point_index].eigenvectors()(row, col);
	}	
	
	double abs2eigenvector(int mesh_point_index, int row, int col) {
		assert(m_eig.size() != 0 && "Run a diagonalization first!");
		return m_eig[mesh_point_index].eigenvectors().cwiseAbs2()(row, col);
	}	

	void push_back(saes_t& saes) {
		m_eig.push_back(saes);	
	}

private:
	std::vector<saes_t> m_eig;
	
};










#endif /* Eigmesh_hpp */
