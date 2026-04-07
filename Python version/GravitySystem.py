import numpy as np

earth_mass = 5.97e24
solar_mass = 1.99e39
"""
The gravitational constant is scaled so that in a system
with a central force potential with the strength of the sun, a
body with unit distance from the origin moving with an angular
velocity of 2pi will complete a full circular orbit 
in unit time.
"""
grav_constant = 4*(np.pi**2)/solar_mass

class Body:
    def __init__(self, r, v, m):
        """
        Class holding the data for a body being affected by gravity

        Parameters
        ----------
        r : ndarray 1D
            Position vector of body
        v : ndarray 1D
            Velocity vector of body
        m : double:
            Mass of body
        """
        self.position = r
        self.velocity = v    
        self.mass = m*earth_mass

class System:
    def __init__(self, central_force = True, interactions = False, central_mass = 1):
        """
        The system that we are simulating

        Parameters
        ----------
        central_force : bool
            Boolean that determines if the system has a central force potential
            located at [0.0, 0.0].
        interactions : bool
            Determines wether the bodies in the system interact with each other via
            gravity or not.
        central_mass : double
            The mass of the body in the center if the central force is turned on,
            affects the strength of the potential. In units of solar mass.
        """
        self.bodies = np.array([])
        self.central_mass = central_mass*solar_mass
        self.central_force_switch = central_force
        self.interactions_switch = interactions

    def add_body(self, body):
        """
        Adds a body to the system

        Parameters
        ----------
        body : Body
            The body that is being added
        """
        self.bodies = np.append(self.bodies, body)

    def central_force(self, i):
        """
        The force that acts upon all bodies if the central force is turned on

        Parameters
        ----------
        i : int
            Picks out the body from the bodies array that is being acted upon by the
            central force.
        
        Returns
        -------
        central_force : ndarray 1D
            Vector of the force that acts upon a body due to the central force potential.
        """
        distance = np.linalg.norm(self.bodies[i].position)
        mass = self.bodies[i].mass
        position = self.bodies[i].position
        central_force = -grav_constant*self.central_mass*mass*position/(distance**3)
        return central_force
    
    def force_body(self, i, j):
        """
        The force that the bodies act upon one another

        Parameters
        ----------
        i : int
            Picks out the body that will have a force acted on it
        j : int
            Picks out the body that is exerting its gravitational
            force on body i.

        Returns
        -------
        force_body : ndarray 1D
            The force vector due to interactions between bodies
        """
        distance = np.linalg.norm(self.bodies[i].position - self.bodies[j].position)
        mass_i = self.bodies[i].mass
        mass_j = self.bodies[j].mass
        position_i = self.bodies[i].position
        position_j = self.bodies[j].position
        force_body = -grav_constant*mass_i*mass_j*(position_i - position_j)/(distance**3)
        return force_body
    
    def total_force_body(self, i):
        """
        The total force that is exerted on a body from all other bodies
        in the system due to gravity.

        Parameters
        ----------
        i : int
            Picks out the body that has the forces acted on it.

        Returns
        -------
        total_force_body : ndarray 1D
            The total force due to gravitational interaction from the
            other bodies.
        """
        total_force_body = np.array([0.0, 0.0])
        for j in range(len(self.bodies)):
            if i != j:
                total_force_body = total_force_body + self.force_body(i, j)

        return total_force_body
    
    def total_force(self, i):
        """
        Returns the total force acting on a body by checking wether or not
        the central potential and interactions are turned on or not.

        Parameters
        ----------
        i : int
            Picks out the body that has the force acted on it.
        
        Returns
        -------
        total_force : ndarray 1D
            The total force that acts on a body
        """
        total_force = np.array([0, 0])
        if self.central_force_switch:
            if self.interactions_switch:
                total_force = self.central_force(i) + self.total_force_body(i)
            else:
                total_force = self.central_force(i)
        elif self.interactions_switch:
            total_force = self.total_force_body(i)
        
        return total_force
    
    def evolve_forward_RK4(self, dt):
        """
        4th order Runge-Kutta method for solving the 2nd order ODE.
        Updates the position vector for each body of the system.

        Parameters
        ----------
        dt : float
            The time step used to solve the equation.
        """
        for i in range(len(self.bodies)):
            #
            temporary_body = Body(self.bodies[i].position, self.bodies[i].velocity, self.bodies[i].mass/earth_mass)
            m_i = self.bodies[i].mass
            
            k_r1 = self.bodies[i].velocity
            k_v1 = self.total_force(i)/m_i
            self.bodies[i].position = self.bodies[i].position + 0.5*dt*k_r1
            self.bodies[i].velocity = self.bodies[i].velocity + 0.5*dt*k_v1
            
            k_r2 = self.bodies[i].velocity
            k_v2 = self.total_force(i)/m_i
            self.bodies[i].position = self.bodies[i].position + 0.5*dt*k_r2
            self.bodies[i].velocity = self.bodies[i].velocity + 0.5*dt*k_v2
            
            k_r3 = self.bodies[i].velocity
            k_v3 = self.total_force(i)/m_i
            self.bodies[i].position = self.bodies[i].position + dt*k_r3
            self.bodies[i].velocity = self.bodies[i].velocity + dt*k_v3
            
            k_r4 = self.bodies[i].velocity
            k_v4 = self.total_force(i)/m_i
            self.bodies[i].position = temporary_body.position + dt*(k_r1 + 2*k_r2 + 2*k_r3 + k_r4)/6
            self.bodies[i].velocity = temporary_body.velocity + dt*(k_v1 + 2*k_v2 + 2*k_v3 + k_v4)/6
    def evolve_forward_RK4_time(self, time_span, dt):
        """
        Evolves the system forward a given amount of time

        Parameters
        ----------
        time_span : float
            The time span that the system evolves in
        dt : float
            The time step taken for evolving the system
        
        Returns
        -------
        """
        counts = int(time_span/dt)
        time_vector = np.zeros(counts)
        position_vector = np.zeros([len(self.bodies), counts, 2])
        body_vector = np.zeros([counts, len(self.bodies)])
        for i in range(len(self.bodies)):
            body_vector[0, i] = self.bodies[i]
        #body_vector[0] = self.bodies[0]
        time_point = 0 #give different name
        for i in range(counts):
            time_vector[i] = time_point
            self.evolve_forward_RK4(dt)
            time_point += dt
        return position_vector