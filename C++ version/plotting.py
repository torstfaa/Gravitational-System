import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

data_collected = np.loadtxt("output.txt")

plt.plot(data_collected[:, 0], data_collected[:, 1])
plt.show()

"""fig = plt.figure(figsize=(8, 6))
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

animation = FuncAnimation(
                    fig=fig,
                    func=update_data,
                    frames=len(position_list[0, 0]),
                    interval=1,
                    repeat=False
)

plt.grid()
animation.save("animation.gif", fps=30)
plt.show()"""