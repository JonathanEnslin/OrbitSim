import matplotlib.pyplot as plt
import matplotlib.animation as animation

from time import time
import numpy as np
from scipy import constants

from simulator import *
from standard_objects import *


m1 = 0.25*2e30
m2 = 2e30

R = 1e11
centre_of_mass = (m1 * 0 + m2 * R) / (m1 + m2)
r1 = centre_of_mass
r2 = R - centre_of_mass

Fg = (constants.gravitational_constant * m1 * m2) / R**2

v1 = np.sqrt(Fg * r1 / m1)
v2 = np.sqrt(Fg * r2 / m2)

star1 = AstronomicalObject("Star 1",
                            mass=m1,
                            position=np.array([0., r1]),
                            velocity=np.array([-v1, 0.]),
                            radius=sun.radius/2,
                            colour = 'y',
)

star2 = AstronomicalObject("Star 2",
                            mass=m2,
                            position=np.array([0., -r2]),
                            velocity=np.array([v2, 0.]),
                            radius=sun.radius,
                            colour = 'y',
)


# Setup of simulation
sim = SimulationSystem()
sun_clone = sun.clone()
earth_clone = earth.clone()
moon_clone = moon.clone()

sim.add_object(sun_clone)
sim.add_object(earth_clone)
sim.add_object(moon_clone)

# sim.add_object(star1.clone())
# sim.add_object(star2.clone())

start_time = time()

# Run simulation
# positions, velocities, time_steps = sim.simulate(dt=step_size, num_steps = int(timespan / step_size), dimensions=2, print_progress=True, track=False)

timespan = 1*constants.year
step_size = 360
print("Step size:", get_time_str(step_size))
positions, velocities, time_steps = sim.simulate(dt=step_size, num_steps = int(timespan / step_size), dimensions=2, print_progress=True)
sim.step(step_size)

end_time = time()
print("Time elapsed: ", end_time - start_time, "seconds")

if input("Would you like the plot to be animated? (y/n) ") == 'y':
    file_name = input("Enter file name, extension will automatically be .mp4: ")
# if True:
    # ani = plot_animated_2D(positions, velocities, time_steps, time_window= 60 * constants.day, frames = 400, interval=15, size_scaling='fixed', plotbounds=((-0.2*constants.au, 1.2*constants.au), (-0.01e12, 1.1e12)), blit=True, name_anno_fontsize=6)
    # ani = plot_animated_2D(positions, velocities, time_steps, time_window= 1*constants.year, frames = 800, interval=15, size_scaling='fixed', blit=True, name_anno_fontsize=6)
    ani = plot_animated_2D(positions, velocities, time_steps, time_window=20*constants.week, frames = 720, interval=15, size_scaling='auto', blit=True, margin_size=0.1, upper_plotbounds=(None, (None, 3e11)), lower_plotbounds=(None, (-0.1e11, None)), time_anno_fontsize=6, name_anno_fontsize=6, ignore_radius=False)
    ani.save(file_name + ".mp4", writer='ffmpeg', fps=24, dpi=250, bitrate=2500, progress_callback=progress_callback_func_factory())
    print()
    # print("Saved")
    # ani = plot_animated_2D(positions, velocities, time_steps, time_window=20*constants.week, frames = 720, interval=15, size_scaling='fixed', blit=True, margin_size=0.1, upper_plotbounds=(None, (None, 3e11)), lower_plotbounds=(None, (-0.1e11, None)), time_anno_fontsize=6, name_anno_fontsize=6, ignore_radius=False)
    # ani.save('moon-around-earth-sun-removed-fixed.mp4', writer='ffmpeg', fps=24, dpi=250, bitrate=2500, progress_callback=progress_callback_func_factory())
    # print()
    # print("Saved")
    # ani = plot_animated_2D(positions, velocities, time_steps, time_window=0.5*constants.week, frames = 720, interval=15, size_scaling='auto', blit=True, margin_size=0.1, upper_plotbounds=(None, (None, 3e11)), lower_plotbounds=(None, (-0.1e11, None)), time_anno_fontsize=6, name_anno_fontsize=6, ignore_radius=False)
    # ani.save('moon-around-earth-sun-removed-half-week window.mp4', writer='ffmpeg', fps=24, dpi=250, bitrate=2500, progress_callback=progress_callback_func_factory())
    # print()
    # print("Saved")

else:
    plot_sim_positions(positions, ignore_radius=True)
    plt.show()
