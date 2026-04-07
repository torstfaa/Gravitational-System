from GravitySystem import *
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

solarSystem = System()
solarSystem.add_body(Body(np.array([1.0, 0.0]), np.array([0.0, 2*np.pi]), 1))
solarSystem.add_body(Body(np.array([2.0, 0.0]), np.array([0.0, np.sqrt(2)*np.pi]), 1))

dt = 0.001
t_end = 1
t_list = np.arange(0, 1, dt)

position_list = np.zeros([solarSystem.bodies.size, 2, t_end/dt])
#for i in range(len(solarSystem.bodies())):


for i in range(len(t_list)):
    solarSystem.evolve_forward_RK4(dt)
    for j in range(len(solarSystem.bodies)):
        position_list[j, :, i] = solarSystem.bodies[j].position
        position_list[j, :, i] = solarSystem.bodies[j].position

interactionSystem = System(central_force = False, interactions = True)
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
axis.scatter([0, 0], [0, 0], color='y') #The 'sun' in the origin
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
plt.show()