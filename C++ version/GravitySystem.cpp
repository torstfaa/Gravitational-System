#include "test_header.h"

vector<double> operator+(const vector<double>& v1, vector<double> v2) {
	vector<double> v3(v1.size());
	for (int i = 0; i < v1.size(); i++) {
		v3[i] = v1[i] + v2[i];
	}

	return v3;
};

vector<double> operator-(vector<double> v1, vector<double> v2) {
	vector<double> v3(v1.size());
	for (int i = 0; i < v1.size(); i++) {
		v3[i] = v1[i] - v2[i];
	}

	return v3;
};

vector<double> operator*(double c, vector<double> v) {
	vector<double> u(v.size());
	for (int i = 0; i < v.size(); i++) {
		u[i] = c*v[i];
	}
	return u;
}

vector<double> operator/(vector<double> v, double c) {
	vector<double> u(v.size());
	for (int i = 0; i < v.size(); i++) {
		u[i] = v[i]/c;
	}
	return u;
}

double operator*(vector<double> v1, vector<double> v2) {
    double return_value = 0;
    for (int i = 0; i < v1.size(); i++) {
        return_value += v1[i]*v2[i];
    }
    return return_value;
}


Body::Body(vector<double> r, vector<double> v, double m) {
    Body::position = r;
    Body::velocity = v;
    Body::mass = m*earth_mass;
}

System::System(bool central_force, bool interactions, double central_mass) {
    System::central_force_switch = central_force;
    System::interactions_switch = interactions;
    System::central_mass = central_mass*solar_mass;
}

void System::add_body(Body body) {
    bodies.push_back(body);
}

vector<Body> System::get_bodies() {
    return bodies;
}

vector<double> System::central_force(int i) {
    double mass = bodies[i].mass;
    vector<double> position = bodies[i].position; //rename?
    double distance = hypot(position[0], position[1]); //Possibly temporary solution
    double force_x = -grav_constant*central_mass*mass*position[0]/(pow(distance, 3));   
    double force_y = -grav_constant*central_mass*mass*position[1]/(pow(distance, 3)); 

    vector<double> force = -grav_constant*central_mass*mass*position/(pow(distance, 3));

    vector<double> central_force = {force_x, force_y};
    return central_force;
}

vector<double> System::force_body(int i, int j) {
    vector<double> position_i = bodies[i].position;
    vector<double> position_j = bodies[j].position;
    double distance = hypot(bodies[i].position[0] - bodies[j].position[0], bodies[i].position[1] - bodies[j].position[1]);
    double mass_i = bodies[i].mass;
    double mass_j = bodies[j].mass;
    double force_x = -grav_constant*mass_i*mass_j*(bodies[i].position[0] - bodies[j].position[0])/(pow(distance, 3));
    double force_y = -grav_constant*mass_i*mass_j*(bodies[i].position[1] - bodies[j].position[1])/(pow(distance, 3));
    vector<double> force_body = {force_x, force_y};
    return force_body;
}

vector<double> System::total_force_body(int i) {
    vector<double> total_force = {0, 0};
    for (int j = 0; j < bodies.size(); j++) {
        if (i != j) {
            total_force = total_force + force_body(i, j);
        }
    }
    return total_force;
}

vector<double> System::total_force(int i) {
    vector<double> total_force = {0.0, 0.0};
    if (central_force_switch) {
        if (interactions_switch) {
            return central_force(i) + total_force_body(i);
        } else {
            return central_force(i);
        }
    } else if (interactions_switch) {
        return total_force_body(i);
    } else {
        return vector<double>(2, 0.0);
    }
}

void System::evolve_forward_RK4(double dt) {
    for (int i = 0; i < bodies.size(); i++) {
        Body temporary_body = Body(bodies[i].position, bodies[i].velocity, bodies[i].mass/earth_mass);
        double m_i = bodies[i].mass;

        vector<double> k_r1 = bodies[i].velocity;
        vector<double> k_v1 = central_force(i)/m_i;
        
        bodies[i].position = bodies[i].position + 0.5*dt*k_r1;
        bodies[i].velocity = bodies[i].velocity + 0.5*dt*k_v1;
        

        vector<double> k_r2 = bodies[i].velocity;
        vector<double> k_v2 = central_force(i)/m_i;
        bodies[i].position = bodies[i].position + 0.5*dt*k_r2;
        bodies[i].velocity = bodies[i].velocity + 0.5*dt*k_v2;

        vector<double> k_r3 = bodies[i].velocity;
        vector<double> k_v3 = central_force(i)/m_i;
        bodies[i].position = bodies[i].position + dt*k_r3;
        bodies[i].velocity = bodies[i].velocity + dt*k_v3;

        vector<double> k_r4 = bodies[i].velocity;
        vector<double> k_v4 = central_force(i)/m_i;
        bodies[i].position = temporary_body.position + dt*(k_r1 + 2*k_r2 + 2*k_r3 + k_r4)/6;
        
        bodies[i].velocity = temporary_body.velocity + dt*(k_v1 + 2*k_v2 + 2*k_v3 + k_v4)/6;
        vector<double> test_vector = temporary_body.position + dt*(k_r1 + 2*k_r2 + 2*k_r3 + k_r4)/6;
    }
}