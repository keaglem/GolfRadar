import numpy
import scipy
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import pytest

# In this package, x is right (aft looking forward), y is up, and z is out
#


gravity_constant = 9.80665


def calculate_air_resistance(drag_coefficient, air_density, surface_area, velocity) :
    """
    Calculating the air resistance on an object

    :param drag_coefficient:
    :param air_density:
    :param surface_area:
    :param velocity:
    :return: air_resistance
    """

    air_resistance = .5 * drag_coefficient * air_density * surface_area * numpy.squeeze(velocity) ** 2

    return air_resistance


def calculate_magnus_force(magnus_coefficient, angular_velocity, air_density, surface_area, velocity):
    """
    Calculating the magnus force on an object

    :param magnus_coefficient:
    :param angular_velocity:
    :param velocity:
    :return: magnus_force
    """

    backspin = velocity.copy()
    #backspin[0] = 0

    sidespin = velocity.copy()
    sidespin[2]=0

    magnus_force = .5 * magnus_coefficient * air_density * surface_area * numpy.cross(numpy.array([1,0,0]), backspin.squeeze() ) ** 2
    magnus_force += .5 * magnus_coefficient * air_density * surface_area * numpy.cross(numpy.array([0,-1,0]), sidespin.squeeze() ) ** 2

    return magnus_force


def calculate_gravity(mass):
    """

    :param mass:
    :return:
    """

    return mass * gravity_constant * numpy.array([0, -1, 0])


def calculate_air_density(temperature_C, relative_humidity=0, pressure_Pa=101325.0):
    """

    :param temperature_C:
    :param relative_humidity:
    :param pressure_Pa:
    :return:
    """

    specific_gas_constant_air = 287.058 # J/ (kg * K)
    specific_gas_constant_water = 461.495 # J/ (kg * K)
    temperature_K = temperature_C + 237.15

    saturation_vapor_pressure_water = 6.1078 * 10 ** (7.5 * temperature_C / temperature_K) / 100

    vapor_pressure_water = relative_humidity * saturation_vapor_pressure_water

    density_air = pressure_Pa / (specific_gas_constant_air * temperature_K)
    density_water = vapor_pressure_water / (specific_gas_constant_water * temperature_K)

    density_humid_air = density_air + density_water

    return density_humid_air


def calculate_reynolds_number(velocity, object_length, air_density):
    """

    :param velocity:
    :param object_length:
    :param air_density:
    :return:
    """

    dynamic_viscosity_air = 1.983e-5

    reynolds_number = air_density * velocity * object_length / dynamic_viscosity_air

    return reynolds_number


def test_reynolds_number():

    test_velocity = numpy.array([[34, 0, 3]]).T  # meters/sec
    test_length = .01  # meters
    test_temperature = 21.  # celsius
    test_relative_humidity = .25  # fraction (0-1)
    test_pressure = 90097.16  # Pa
    test_air_density = calculate_air_density(test_temperature,
                                             test_relative_humidity,
                                             test_pressure)
    reynolds_number = calculate_reynolds_number(test_velocity,
                                                test_length,
                                                test_air_density)
    print('Reynolds number: {}'.format(reynolds_number))


def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = numpy.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = numpy.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = numpy.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([0, 2*plot_radius])


def test_ball_flight():

    test_drag_coefficient = .2  #
    test_magnus_coefficient = .17
    test_exit_speed = 67  # m/s
    test_spin = numpy.array([[0,
                              0,
                              -150]]).T  # rad/s
    test_exit_launch_angle = 20 * numpy.pi / 180.0
    test_velocity = numpy.array([[0,
                                  test_exit_speed*numpy.sin(test_exit_launch_angle),
                                  test_exit_speed*numpy.cos(test_exit_launch_angle)]]).T  # meters/sec
    test_length = .042  # meters
    surface_area_ball = 4 * numpy.pi * (test_length / 2) ** 2
    test_temperature = 21.  # celsius
    test_relative_humidity = .25  # fraction (0-1)
    test_pressure = 90097.16  # Pa
    ball_mass = .0433  # kg
    initial_position = numpy.array([[0, 0, 0]]).T  # initial position in meters

    time_vector = numpy.linspace(start=0, stop=30, num=50000)

    current_position = numpy.zeros((3,len(time_vector)))
    current_velocity = numpy.zeros((3,len(time_vector)))

    current_velocity[:,0] = test_velocity[:,0]
    current_position[:,0] = initial_position[:,0]

    delta_time = time_vector[1] - time_vector[0]

    for index, current_time in enumerate(time_vector[:-1]):

        test_air_density = calculate_air_density(test_temperature,
                                                 test_relative_humidity,
                                                 test_pressure)
        resistance = -calculate_air_resistance(test_drag_coefficient,
                                              test_air_density,
                                              surface_area_ball,
                                              current_velocity[:,index])
        magnus_force = calculate_magnus_force(test_magnus_coefficient,
                                              test_spin,
                                              test_air_density,
                                              surface_area_ball,
                                              current_velocity[:,index])

        resistance += magnus_force
        resistance += calculate_gravity(ball_mass)
        resistance /= ball_mass

        accel_delta_time = (time_vector[index+1] ** 2 - time_vector[index] ** 2)

        current_velocity[:,index+1] = current_velocity[:,index] + \
                                    resistance * delta_time

        current_position[:,index+1] += current_velocity[:,index] * delta_time + \
                                     current_position[:,index]



    filtered_ball_position = current_position[:, current_position[1,:] > 0]
    fig = pyplot.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot(filtered_ball_position[0, :],
            filtered_ball_position[2, :],
            filtered_ball_position[1, :])

    pyplot.axis('equal')
    pyplot.xlabel('Left/Right')
    pyplot.ylabel('In/Out')
    set_axes_equal(ax)
    ax.set_zlabel('Up/Down')
    ax.plot(filtered_ball_position[0, :],
            filtered_ball_position[2, :], 'g--')

    pyplot.show()