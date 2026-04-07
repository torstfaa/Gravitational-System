#include<vector>
#include<cmath>
#include<numeric>
#include<iostream>

using namespace std;

constexpr double pi = 3.14159265;
constexpr double earth_mass = 5.97e24;
constexpr double solar_mass = 1.99e39;
const double grav_constant = 4*pow(pi, 2)/solar_mass;

vector<double> operator+(const vector<double>& v1, vector<double> v2);
vector<double> operator-(vector<double> v1, vector<double> v2);
vector<double> operator*(double c, vector<double> v);
vector<double> operator/(vector<double> v, double c);

struct Body {
    vector<double> position;
    vector<double> velocity;
    double mass;
    Body(vector<double> r, vector<double> v, double m);
};

class System {
    private:
        vector<Body> bodies;
        bool central_force_switch;
        bool interactions_switch;
        double central_mass;
    public:
        void add_body(Body body);
        vector<Body> get_bodies();
        vector<double> central_force(int i);
        vector<double> force_body(int i, int j);
        vector<double> total_force_body(int i);
        vector<double> total_force(int i);
        void evolve_forward_RK4(double dt);
        System(bool central_force = true, bool interactions = false, double central_mass = 1);
};

