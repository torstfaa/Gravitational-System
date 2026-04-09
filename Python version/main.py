from GravitySystem import *
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import time

"""
Here are provided multiple examples for how certain
systems evolve.

Each example saves as a gif which are already provided.
"""

"""
Example 1:
Two non-interacting bodies moving in circular orbits in
a central force potential.
"""

solarSystem = System()
solarSystem.add_body(Body(np.array([1.0, 0.0]), np.array([0.0, 2*np.pi]), 1))
solarSystem.add_body(Body(np.array([0.0, 2.0]), np.array([-np.sqrt(2)*np.pi, 0]), 1))

dt = 0.0001
t_end = 4
steps = int(t_end/dt)

position_list = np.zeros([solarSystem.bodies.size, 2, steps])
#Update the position list for initial positions: 
for i in range(len(solarSystem.bodies)):
    position_list[i, 0, 0] = solarSystem.bodies[i].position[0]
    position_list[i, 1, 0] = solarSystem.bodies[i].position[1]

start_time = time.time()
for i in range(1, steps):
    solarSystem.evolve_forward_RK4(dt)
    for j in range(len(solarSystem.bodies)):
        position_list[j, :, i] = solarSystem.bodies[j].position
        position_list[j, :, i] = solarSystem.bodies[j].position
print("Simulation time: " + str(time.time() - start_time))

fig = plt.figure(figsize=(8, 8))
axis = fig.add_subplot()
axis.scatter([0, 0], [0, 0], color='y') #The 'sun' in the origin
axis.set_xlim([min(position_list[1, 0]) - 0.5, max(position_list[1, 0]) + 0.5])
axis.set_ylim([min(position_list[1, 1]) - 0.5, max(position_list[1, 1]) + 0.5])
axis.set_title("Two non-interacting bodies moving in circular orbit")

path1, = axis.plot([], [], color='k') #The path of the body
body1, = axis.plot([], [], 'o', markersize=8, color='b') #The body itself
path2, = axis.plot([], [], color='k') #The path of the body
body2, = axis.plot([], [], 'o', markersize=8, color='r') #The body itself
body1.set_label("Body 1")
body2.set_label("Body 2")
axis.legend()

"""
Here we create a new list with only every tenth element.
This is done to give the animation fewer frames, thus
speeding it up. Otherwise the animation, while functional,
takes a bit too long to play out.
"""
new_position_list = position_list[:, :, ::50]
start_time = time.time()

def update_data(frame):
    path1.set_data(new_position_list[0][0][:frame], new_position_list[0][1][:frame])
    body1.set_data([new_position_list[0][0][frame]], [new_position_list[0][1][frame]])
    path2.set_data(new_position_list[1][0][:frame], new_position_list[1][1][:frame])
    body2.set_data([new_position_list[1][0][frame]], [new_position_list[1][1][frame]])

    return path1,

animation = FuncAnimation(fig=fig, func=update_data, frames=len(new_position_list[0, 0]), interval=25, repeat=False)
print("Plotting time: " + str(time.time() - start_time))

plt.grid()
animation.save("solarSystem.gif")
plt.show()

"""
Example 2:
Two interacting bodies without a central force
potential. One has the mass of the sun, the other
one tenth of it. 
"""

"""interactionSystem = System(central_force = False, interactions = True)
interactionSystem.add_body(Body(np.array([1.0, 0.0]), np.array([0.0, 1]), 0.1*solar_mass/earth_mass))
interactionSystem.add_body(Body(np.array([-1.0, 0.0]), np.array([0.0, -1]), solar_mass/earth_mass))

dt = 0.001
t_list = np.arange(0, 1, dt)

position_list = np.zeros([interactionSystem.bodies.size, 2, t_list.size])
velocity_list = np.zeros_like(position_list)
for i in range(len(position_list)):
    position_list[i, :, 0] = interactionSystem.bodies[i].position

for i in range(len(t_list)):
    interactionSystem.evolve_forward_RK4(dt)
    for j in range(len(interactionSystem.bodies)):
        position_list[j, :, i] = interactionSystem.bodies[j].position
        position_list[j, :, i] = interactionSystem.bodies[j].position

fig = plt.figure(figsize=(8, 6))
axis = fig.add_subplot()
axis.set_xlim([min(position_list[0, 0]) - 0.5, max(position_list[0, 0]) + 0.5])
axis.set_ylim([min(position_list[0, 1]) - 0.5, max(position_list[0, 1]) + 0.5])
axis.set_title("A dummy title")

animated_plot, = axis.plot([], [], color='k')
animated_plot2, = axis.plot([], [], 'o', markersize=8, color='b')
animated_plot3, = axis.plot([], [], color='k')
animated_plot4, = axis.plot([], [], 'o', markersize=8, color='r')

def update_data(frame):
    animated_plot.set_data(position_list[0][0][:frame], position_list[0][1][:frame])
    animated_plot2.set_data([position_list[0][0][frame]], [position_list[0][1][frame]])
    animated_plot3.set_data(position_list[1][0][:frame], position_list[1][1][:frame])
    animated_plot4.set_data([position_list[1][0][frame]], [position_list[1][1][frame]])

    return animated_plot,

animation = FuncAnimation(fig=fig, func=update_data, frames=len(position_list[0, 0]), interval=25, repeat=False)

plt.grid()
animation.save("animation.gif")
plt.show()"""

"""
Example 3:
Scattering of multiple bodies on a central
force potential.
"""