#include<fstream>
#include "test_header.h"
//#include <armadillo>

int main() {
	vector<double> r1 = {1, 0};
	vector<double> v1 = {0, 2*pi};
	Body b1(r1, v1, 1);

	System s1(true, false, 1);
	s1.add_body(b1);
	double dt = 0.0005;
	double time_span = 1;
	double time_point = 0;
	vector<double> t_list(int(time_span/dt));
	t_list[0] = time_point;
	vector<vector<vector<double>>> position_list(s1.get_bodies().size(), vector<vector<double>>(t_list.size(), vector<double>(2)));
	position_list[0][0][0] = 1;
	position_list[0][0][1] = 0;

	for (int i = 1; i < t_list.size(); i++) {
		time_point += dt;
		t_list[i] = time_point;
		s1.evolve_forward_RK4(dt);
		for (int j = 0; j < s1.get_bodies().size(); j++) {
			position_list[j][i][0] = s1.get_bodies()[j].position[0];
			position_list[j][i][1] = s1.get_bodies()[j].position[1];
		}
	}

	ofstream myFile("output.txt");

	for (int i; i < position_list[0].size(); i++) {
		myFile << position_list[0][i][0] << '\t' << position_list[0][i][1] << '\n';
	}

	return 0;
}

