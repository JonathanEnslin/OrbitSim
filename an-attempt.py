import numpy as np
from scipy import constants
from time import sleep, time
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def circular_orbital_velocity(mass, orbital_radius):
    return np.sqrt(constants.gravitational_constant * mass / orbital_radius)

def orbital_radius(mass, velocity):
    return constants.gravitational_constant * mass / velocity**2

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or np.fabs(value - array[idx-1]) < np.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

def find_nearest_i(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or np.fabs(value - array[idx-1]) < np.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

def get_time_str(curr_time):
    if curr_time > constants.year:
        time_str = f"T: {round(curr_time / constants.year, 3)} years"
    elif curr_time > constants.week:
        time_str = f"T: {round(curr_time / constants.week, 3)} weeks"
    elif curr_time > constants.day:
        time_str = f"T: {round(curr_time / constants.day, 3)} days"
    elif curr_time > constants.hour:
        time_str = f"T: {round(curr_time / constants.hour, 3)} hours"
    elif curr_time > constants.minute:
        time_str = f"T: {round(curr_time / constants.minute, 3)} minutes"
    else:
        time_str = f"T: {round(curr_time, 3)} seconds"
    return time_str


def get_v_str(v):
    if abs(v) > 1000:
        v_str = f"{round(v / 1000, 2)} km/s"
    else:
        v_str = f"{round(v, 2)} m/s"
    return v_str

def gravitational_force(m1, m2, r):
    return constants.gravitational_constant * m1 * m2 / r**2

def euclidean_distance(point1, point2):
    return np.linalg.norm(point1 - point2)

class AstronomicalObject:
    def __init__(self, name, mass, position, velocity=np.array([0., 0.]), radius=None, colour='b', colour2='black'):
        self.name = name
        self.mass = mass
        self.position = position  # a position matrix
        self.velocity = velocity  # a velocity vector
        self.radius = radius
        self.colour = colour
        self.colour2 = colour2

    def __str__(self):
        return f'Name: {self.name}\nMass: {self.mass}kg\nPosition: {self.position}m\nVelocity: <{self.velocity[0]} {self.velocity[1]}>m/s\nRadius: {self.radius/1000}km'

    def step(self, dt, acceleration=np.array([0.0, 0.0]), accel_disp = True):
        # Update position
        self.position = self.position + self.velocity * dt
        # apply displacement due to acceleration
        if accel_disp:
            self.position +=  0.5 * acceleration * dt**2
        # Update velocity
        self.velocity = self.velocity + acceleration * dt
        return self.position, self.velocity

    def movement_str(self):
        return f'Position: {self.position}m  Velocity: <{self.velocity[0]} {self.velocity[1]}>m/s'

    def __hash__(self):
        return hash(self.name)

    def clone(self):
        # shallow copies, will need to change if we ever store objects in the position or velocity
        return AstronomicalObject(self.name, self.mass, np.copy(self.position), np.copy(self.velocity), self.radius, self.colour, self.colour2)

    def create_variant(self, **kwargs):
        new_obj = self.clone()
        new_obj.set_properties(**kwargs)
        new_obj.velocity = new_obj.velocity.copy()
        new_obj.position = new_obj.position.copy()
        return new_obj

    def set_properties(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)


class SimulationSystem:
    def __init__(self):
        self.astro_objects = []
        self.time = 0


    def add_object(self, astronomical_object: AstronomicalObject):
        self.astro_objects.append(astronomical_object)


    def step(self, dt):
        # TODO: @JonathanEnslin, this method is extreeeeemely inneficient, 1. change to not do duplicate calcs, 2. use vectorisation
        # step all objects
        self.time += dt
        if len(self.astro_objects) == 0:
            return
        # Start with 0 resulting force on all objects
        resultant_forces = [np.zeros(self.astro_objects[0].position.shape, dtype=np.float64) for _ in range(len(self.astro_objects))]
        for astro_i, astro_obj in enumerate(self.astro_objects):
            for other_astro_i, other_astro_obj in enumerate(self.astro_objects[astro_i:], astro_i):
                if astro_i == other_astro_i:
                    # Force should be zero
                    resultant_forces[astro_i] += np.zeros(other_astro_obj.position.shape)
                    continue
                # Calculate displacement matrix
                displacement_matrix = other_astro_obj.position - astro_obj.position
                # Calculate euclidean distance between two
                distance = np.linalg.norm(displacement_matrix)
                # print("Distance: ", distance)
                # Calcultate the component ratios
                component_ratios = displacement_matrix / distance
                # Calculate gravitational pull between the two
                Fg = gravitational_force(astro_obj.mass, other_astro_obj.mass, distance)
                # Get force components for astro_obj
                Fg_on_astro_obj = Fg * component_ratios
                # Add resultant forces
                resultant_forces[astro_i] += Fg_on_astro_obj
                resultant_forces[other_astro_i] -= Fg_on_astro_obj

        # Apply steps
        for astro_obj, resultant_force in zip(self.astro_objects, resultant_forces):
            # print("Resultant Force on:", astro_obj.name, resultant_force)
            # print("Acceleration on:   ", astro_obj.name, resultant_force / astro_obj.mass)

            astro_obj.step(dt, resultant_force / astro_obj.mass)


    def simulate(self, dt, num_steps, dimensions=2, print_progress = False):
        # TODO: @JonathanEnslin, use dicts instead and map using object name, additionally, look at include time steps
        # Initialise position and velocity tracking arrays
        positions = {astro_obj: [astro_obj.position] for astro_obj in self.astro_objects}
        velocities = {astro_obj: [astro_obj.velocity] for astro_obj in self.astro_objects}
        # ======== Print progress code ========
        curr_time, prev_time = time(), time()
        total_elapsed_time_between_prints = 1 # init to 1 to avoid zero division
        num_prints = 0
        print_period = 100 # Just a random number
        # =====================================
        time_steps = [self.time]
        for i in range(num_steps):
            if print_progress and i % print_period == 0:
                num_prints += 1
                curr_time = time()
                time_elapsed_between_prints = curr_time - prev_time
                total_elapsed_time_between_prints += time_elapsed_between_prints
                average_elapsed_time_betwen_prints = total_elapsed_time_between_prints / num_prints
                periods_left = (num_steps - i) / print_period
                estimated_time_left = periods_left * average_elapsed_time_betwen_prints
                prev_time = curr_time
                print(f"Progress: {round(i / num_steps * 100, 2)}%  -  ETA: {get_time_str(estimated_time_left)}        ", end='\r')
            self.step(dt)
            # Add updated positions and velocities
            for astro_obj in self.astro_objects:
                positions[astro_obj].append(astro_obj.position)
                velocities[astro_obj].append(astro_obj.velocity)
            time_steps.append(self.time)

        # limit the velocity and position to dimensions dimensions
        for astro_obj in self.astro_objects:
            positions[astro_obj] = [pos[:dimensions] for pos in positions[astro_obj]]
            velocities[astro_obj] = [vel[:dimensions] for vel in velocities[astro_obj]]

        if print_progress:
            print(' '*50, "\rProgress: 100%")
            print("Complete.")
        return positions, velocities, time_steps


# Axes to plot should be a list of lists of length 2, where each list contains the axes to plot, e.g. [['x', 'y'], ['x', 'z']]
def plot_sim_positions(positions, fig=None, ax=None, xlim=None, ylim=None, **kwargs):
    fig_or_ax_None = fig is None or ax is None
    ignore_radius = kwargs.get("ignore_radius", False)
    plt.rcParams['figure.figsize'] = [8, 8]


    # ================== plotting ==================
    if fig is None or ax is None:
        fig, (ax) = plt.subplots()

    ax.set_title(f"Solar System")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    data_len = len(positions[list(positions.keys())[0]])
    skip_size = max(int(data_len / 262801), 1)

    for astro, positions in positions.items():
        movement = list(zip(*positions))
        x_mov = movement[0][::skip_size]
        y_mov = movement[1][::skip_size]


        radius = astro.radius/1e5 if astro.radius is not None else None

        ax.plot(x_mov, y_mov, color=astro.colour, linewidth=0.4, alpha=1)
        ax.scatter(*positions[0], c=astro.colour, s=radius if not ignore_radius else None, alpha=0.6)
        ax.scatter(*positions[-1], c=astro.colour2, s=radius if not ignore_radius else None, alpha=0.4)
        ax.annotate(f"{astro.name} initial", positions[0], xytext=(10, 0), textcoords="offset pixels")
        ax.annotate(f"{astro.name} final", positions[-1], xytext=(10, -10), textcoords="offset pixels")

    if fig_or_ax_None:
        ax.set(aspect='equal')
        # ax = plt.gca()
        # ax.set_aspect('equal', adjustable='box')

    plt.show()

# if named par size_scaling is set to auto, it will scale automatically, withing the potential upper and lower bounds passed in in upper_plotbounds -> ((x_min, x_max), (y_min, y_max)) and lower_plotbounds -> ((x_min, x_max), (y_min, y_max))
# if named par size_scaling is set to fixed, it will remain to the fixed size in plotbounds -> ((x_min, x_max), (y_min, y_max))
def plot_animated_2D(positions, velocities, time_steps, time_window=None, frames=500, interval=40, **kwargs):
    velocities = velocities or {}
    blit =  kwargs['blit'] if 'blit' in kwargs else True
    time_anno_fontsize = kwargs['time_anno_fontsize'] if 'time_anno_fontsize' in kwargs else 8
    name_anno_fontsize = kwargs['name_anno_fontsize'] if 'name_anno_fontsize' in kwargs else 8
    size_scaling = kwargs['size_scaling'] if 'size_scaling' in kwargs else 'auto' # Default is 'auto', can be either 'auto' or 'fixed'
    # plotbounds = kwargs['plotbounds'] if 'plotbounds' in kwargs else (None, None) # Default is None, can be either None or ((x_min, x_max), (y_min, y_max)), any value can be None in the list/tuple
    plotbounds = kwargs['plotbounds'] if 'plotbounds' in kwargs else ((None, None), (None, None)) # Default is None, can be either None or ((x_min, x_max), (y_min, y_max)), any value can be None in the list/tuple
    upper_plotbounds = kwargs['upper_plotbounds'] if 'upper_plotbounds' in kwargs else ((None, None), (None, None))
    # upper_plotbounds = kwargs['upper_plotbounds'] if 'upper_plotbounds' in kwargs else (None, None)
    lower_plotbounds = kwargs['lower_plotbounds'] if 'lower_plotbounds' in kwargs else ((None, None), (None, None))
    # lower_plotbounds = kwargs['lower_plotbounds'] if 'lower_plotbounds' in kwargs else (None, None)
    margin_size = kwargs['margin_size'] if 'margin_size' in kwargs else 0.1
    ignore_radius = kwargs.get("ignore_radius", False)
    # Check if the size scaling is valid
    if size_scaling not in ['auto', 'fixed']:
        raise ValueError(f"size_scaling must be either 'auto' or 'fixed', not {size_scaling}")

    for bound, name in zip([plotbounds, upper_plotbounds, lower_plotbounds], ['plotbounds', 'upper_plotbounds', 'lower_plotbounds']):
        if not isinstance(bound, tuple) or len(bound) != 2:
            raise ValueError(f"{name} must be a tuple of length 2, not {bound}")
        for bound_i, bound in enumerate(bound):
            if not isinstance(bound, tuple) and bound is not None:
                raise ValueError(f"{name}[{bound_i}] must be a tuple of length 2, or None, not {bound}")
            if bound is None:
                continue
            for val in bound:
                if not isinstance(val, (int, float)) and val is not None:
                    raise ValueError(f"Bound values in {name} must be ints or floats, not {val} of type {type(val)}")

    # convert plotbounds to list, and inner to tuples to lists as well
    plotbounds = list(plotbounds)
    for i, bound in enumerate(plotbounds):
        if bound is not None:
            plotbounds[i] = list(bound)

    # convert upper_plotbounds to list, and inner to tuples to lists as well
    upper_plotbounds = list(upper_plotbounds)
    for i, bound in enumerate(upper_plotbounds):
        if bound is not None:
            upper_plotbounds[i] = list(bound)

    # convert lower_plotbounds to list, and inner to tuples to lists as well
    lower_plotbounds = list(lower_plotbounds)
    for i, bound in enumerate(lower_plotbounds):
        if bound is not None:
            lower_plotbounds[i] = list(bound)


    velocities_seperated = {astro: list(zip(*velocities[astro])) for astro in velocities}

    # Set flag for autoscaling
    autoscaling = size_scaling == 'auto'

    # Set the plot bounds default values to the max and x and y and min x and y values out of all the lists/arrays in the positions dict
    overall_min_y = min([min([pos[1] for pos in pos_list]) for pos_list in positions.values()])
    overall_max_y = max([max([pos[1] for pos in pos_list]) for pos_list in positions.values()])
    overall_min_x = min([min([pos[0] for pos in pos_list]) for pos_list in positions.values()])
    overall_max_x = max([max([pos[0] for pos in pos_list]) for pos_list in positions.values()])

    if not autoscaling:
        # If the plotbounds are not set, set them to the overall min and max values
        if plotbounds[0] is None:
            plotbounds = [[overall_min_x - margin_size * np.abs(overall_min_x), overall_max_x + margin_size * np.abs(overall_max_x)], plotbounds[1]]
        if plotbounds[1] is None:
            plotbounds = [plotbounds[0], [overall_min_y - margin_size * np.abs(overall_min_y), overall_max_y + margin_size * np.abs(overall_max_y)]]
        if plotbounds[0][0] is None:
            plotbounds[0][0] = overall_min_x - margin_size * np.abs(overall_min_x)
        if plotbounds[0][1] is None:
            plotbounds[0][1] = overall_max_x + margin_size * np.abs(overall_max_x)
        if plotbounds[1][0] is None:
            plotbounds[1][0] = overall_min_y - margin_size * np.abs(overall_min_y)
        if plotbounds[1][1] is None:
            plotbounds[1][1] = overall_max_y + margin_size * np.abs(overall_max_y)
    else:
        # If the plotbounds are not set, set them to the overall min and max values
        # Margins are added if the plotbounds are not explicitly set
        for i, bound in enumerate(upper_plotbounds):
            if bound is None:
                upper_plotbounds[i] = [float('-inf'), float('inf')]

        for i, bound in enumerate(lower_plotbounds):
            if bound is None:
                lower_plotbounds[i] = [-0., -0.]

        for bound in upper_plotbounds:
            if bound[0] is None:
                bound[0] = float('-inf')
            if bound[1] is None:
                bound[1] = float('inf')

        for i in range(len(lower_plotbounds)):
            for j in range(len(lower_plotbounds[i])):
                if lower_plotbounds[i][j] is None:
                    lower_plotbounds[i][j] = -0.

    # Create figure
    fig = plt.figure()
    # axis = plt.axes()
    # axis.set_aspect('equal', 'box')
    axis = plt.axes(xlim=plotbounds[0], ylim=plotbounds[1], aspect='equal')

    # Extract the astronomical objects
    astros = list(positions.keys())
    # Create a line for each astro
    lines = [axis.plot([], [], lw = 1, alpha=0.5, c=astro.colour)[0] for astro in astros]
    # Create a scatter plot for each astro
    scatters = [axis.scatter([], [], c=astro.colour, s=None if (astro.radius is None or ignore_radius) else 2*astro.radius/2e6) for astro in astros]

    # Create the annotations
    annotations = [axis.annotate(f"{astro.name}", (0, 0), alpha=0.7, xytext=(np.random.randint(-10, 10), np.random.randint(-10, 10)), textcoords="offset pixels", fontsize=name_anno_fontsize) for astro in astros]

    # Add speed to the annotations where astro is in seperated_velocities
    for i, astro in enumerate(astros):
        if astro in velocities_seperated:
            annotations[i].set_text(f"{astro.name}\nx_v: {get_v_str(velocities_seperated[astro][0][0])}\ny_v: {get_v_str(velocities_seperated[astro][1][0])}")

    # Create the annotation for time
    anno_time = axis.annotate(f"Time: {time_steps[0]}", xy=(0,0), xycoords='axes fraction', fontsize=time_anno_fontsize,
            xytext=(0,0), textcoords='offset points',
            bbox=dict(boxstyle="square", fc="w"),
            arrowprops=dict(arrowstyle="->"))

    # Extract the x and y positions for all astros
    x_positions = []
    y_positions = []
    for astro in astros:
        x_pos, y_pos = zip(*positions[astro])
        x_positions.append(x_pos)
        y_positions.append(y_pos)

    # get the number of data points
    data_len = len(positions[astros[0]])
#     print(data_len)
    def init():
        return *lines, *scatters

    skip_size = max(int(data_len / 262801), 1)
    # print(skip_size)
    def animate(i):
        xi, yi = int(i*data_len/frames), int(i*data_len/frames)
        if time_window is not None:
            # Only display within time window
            curr_time = time_steps[yi]  # The time step of the current frame
            # Find what time step should be searched for
            search_time = curr_time - time_window
            # Find the index of the time step
            time_i = np.searchsorted(time_steps, search_time, side='right')
        else:
            time_i = 0

        #  set x/y_min and x/y_max to negative and positive infinity
        y_min, y_max = float('inf'), float('-inf')
        x_min, x_max = float('inf'), float('-inf')

        for astro_index in range(len(astros)):
            x_data = x_positions[astro_index][time_i:xi:skip_size]
            y_data = y_positions[astro_index][time_i:yi:skip_size]
            lines[astro_index].set_data(x_data, y_data)
            # scale the plot
            if len(x_data) != 0 and autoscaling:
                temp_y_min = min(y_data)
                temp_y_max = max(y_data)
                temp_x_min = min(x_data)
                temp_x_max = max(x_data)
                # update the min and max
                y_min = min(y_min, temp_y_min)
                y_max = max(y_max, temp_y_max)
                x_min = min(x_min, temp_x_min)
                x_max = max(x_max, temp_x_max)

                # take the upper and lower bounds of the plot into consideration
                if upper_plotbounds is not None: # If statement is not necessary, but is there for clarity/just incase
                    # Not allowed to go bigger than the upper bounds
                    y_max = min(y_max, upper_plotbounds[1][1])
                    x_max = min(x_max, upper_plotbounds[0][1])
                    # Not allowed to go smaller than upper bounds lower limit
                    y_min = max(y_min, upper_plotbounds[1][0])
                    x_min = max(x_min, upper_plotbounds[0][0])
                if lower_plotbounds is not None:
                    # Not allowed to go larger than lower bounds lower limit
                    y_min = min(y_min, lower_plotbounds[1][0])
                    x_min = min(x_min, lower_plotbounds[0][0])
                    # Not allowed to go smaller than lower bounds upper limit
                    y_max = max(y_max, lower_plotbounds[1][1])
                    x_max = max(x_max, lower_plotbounds[0][1])

                # set the figure limits to be equal to the extremes
                # axis.set_ylim(bottom=minmin - 0.5 * np.abs(minmin), top=maxmax + 0.5 * np.abs(maxmax))
                # axis.set_xlim(left=minmin - 0.5 * np.abs(minmin), right=maxmax + 0.5 * np.abs(maxmax))
                axis.set_ylim(bottom=y_min - margin_size * np.abs(y_min), top=y_max + margin_size * np.abs(y_max))
                axis.set_xlim(left=x_min - margin_size * np.abs(x_min), right=x_max + margin_size * np.abs(x_max))

            scatters[astro_index].set_offsets([x_positions[astro_index][xi], y_positions[astro_index][yi]])
            annotations[astro_index].xy = (x_positions[astro_index][xi], y_positions[astro_index][yi])


            if astros[astro_index] in velocities_seperated:
                annotations[astro_index].set_text(f"{astros[astro_index].name}\nx_v: {get_v_str(velocities_seperated[astros[astro_index]][0][xi])}\ny_v: {get_v_str(velocities_seperated[astros[astro_index]][1][xi])}")

        # update the time text
        anno_time.set_text(get_time_str(time_steps[yi]))

        return *lines, *scatters

    return animation.FuncAnimation(fig, animate,
                            init_func = init,
                            frames = frames,
                            interval = interval,
                            blit = blit)


def progress_callback_func_factory():
    start_time = time()
    def progress_callback_func(frame, total_frames):
        frame += 1
        elapsed_time = time() - start_time
        avg_time_per_frame = elapsed_time / frame
        remaining_time = avg_time_per_frame * (total_frames - frame)
        print(f"Frame {frame}/{total_frames} - {(frame/total_frames*100):.2f}% - Remaining time: {get_time_str(remaining_time)} - Elapsed time: {get_time_str(elapsed_time)}", ' '*10, end='\r')
    return progress_callback_func

# Default solar system objects
sun_3D = AstronomicalObject("Sun",
                         mass=1.989e+30,
                         position=np.array([0., 0, 0]),
                         velocity=np.array([0., 0, 0]),
                         radius=696_340_000,
                        #  radius=None,
                         colour = 'y',
                         colour2 = 'orange'
                        )
earth_3D_tilted = AstronomicalObject("Earth Tilted",
                           mass=5.972e+24,
                           position=np.array([np.cos(45*constants.degree)*constants.au, 0.,  np.sin(45*constants.degree)*constants.au]),
                           velocity=np.array([0., 107_226*constants.kmh, 0]),
                           radius=None,
                           colour = 'blue',
                           colour2 = 'blue'
                          )
earth_3D = AstronomicalObject("Earth",
                           mass=5.972e+24*1e4,
                           position=np.array([0., constants.au, 0.]),
                           velocity=np.array([-107_226*constants.kmh, 0., 0]),
                           radius=None,
                           colour = 'green',
                           colour2 = 'green'
                          )
moon = AstronomicalObject("Moon",
                            mass=7.34767309e+22,
                            position=np.array([constants.au + 384_400_000, 0., 0]),
                            velocity=np.array([0., 107_226*constants.kmh + 1_022*constants.kmh, 0]),
                            radius=None,
                            colour = 'grey',
                            colour2 = 'grey'
                            )

sim = SimulationSystem()
sim.add_object(sun_3D)
sim.add_object(earth_3D)
sim.add_object(earth_3D_tilted)

timespan = 20*constants.year
step_size = 1800

positions, velocities, time_steps = sim.simulate(dt=step_size, num_steps = int(timespan / step_size), dimensions=3, print_progress=True)


earth_positions = positions[earth_3D]
earth_positions = list(zip(*earth_positions))

eart_tilted_positions = positions[earth_3D_tilted]
eart_tilted_positions = list(zip(*eart_tilted_positions))

sun_positions = positions[sun_3D]
sun_positions = list(zip(*sun_positions))

# find the max value out of all the values
max_val = max([max(earth_positions[0]), max(earth_positions[1]), max(earth_positions[2])])
max_val = max([max_val, max([max(eart_tilted_positions[0]), max(eart_tilted_positions[1]), max(eart_tilted_positions[2])])])
max_val = max([max_val, max([max(sun_positions[0]), max(sun_positions[1]), max(sun_positions[2])])])
# find the min value out of all the values
min_val = min([min(earth_positions[0]), min(earth_positions[1]), min(earth_positions[2])])
min_val = min([min_val, min([min(eart_tilted_positions[0]), min(eart_tilted_positions[1]), min(eart_tilted_positions[2])])])
min_val = min([min_val, min([min(sun_positions[0]), min(sun_positions[1]), min(sun_positions[2])])])

# get the difference between max and min
diff = max_val - min_val

from mayavi import mlab
import numpy as np




mlab.figure(bgcolor=(0.5, 0.5, 0.5))
mlab.plot3d(*earth_positions, time_steps, colormap='Blues', tube_radius=diff/1000)
mlab.plot3d(*eart_tilted_positions, time_steps, colormap='Greens', tube_radius=diff/1000)
mlab.plot3d(*sun_positions, time_steps, colormap='Reds', tube_radius=diff/1000)
# mlab.plot3d(*earth_positions, color=(0,0,1), tube_radius=diff/1000)
# mlab.plot3d(*eart_tilted_positions, color=(0,1,0), tube_radius=diff/1000)
# mlab.plot3d(*sun_positions, color=(1,0,0), tube_radius=diff/1000)
mlab.show()
