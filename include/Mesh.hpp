///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** Mesh.hpp ***                                 //
//                                                                           //
// Base class responsible for holding discretized mesh points and weights    //
// templatized over dimension of the mesh                                    //
//                                                                           //
// created November 26, 2017                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Mesh_hpp
#define Mesh_hpp

#include <Eigen/Dense>
#include <numeric>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

class Mesh 
{
public:
	Mesh(std::string klist_file) {

		std::ifstream file(klist_file);
		assert(file.is_open());
		
		std::string line;
		Eigen::Vector3d invec;
		double trash, weight, div;

		while (std::getline(file, line)) {
			if (line == "END") break;
			std::stringstream input(line);
			input >> trash;
			for (int i = 0; i < 3; ++i) input >> invec(i);
			input >> div;
			m_points.push_back(invec * 2 * M_PI / div);
			input >> weight;
			m_weights.push_back(weight);
		}	
	}	

	Eigen::Vector3d point(int index) { return m_points[index]; }
	int weight(int index) { return m_weights[index]; }
	int wedge() { return m_points.size(); }
	int full() { return std::accumulate(m_weights.begin(), m_weights.end(), 0); }

private:
	std::vector<Eigen::Vector3d> m_points;
	std::vector<int> m_weights;
};

#endif /* Mesh_hpp */
