import numpy
from scipy import signal
import matplotlib
from GolfBallTracking.physics import physics
from matplotlib import pyplot

# In this package, x is right (aft looking forward), y is up, and z is out


def plot_position_vector(position_vector):
    filtered_ball_position = position_vector[:, position_vector[1, :] > 0]
    fig = pyplot.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot(filtered_ball_position[0, :],
            filtered_ball_position[2, :],
            filtered_ball_position[1, :])

    pyplot.axis('equal')
    pyplot.xlabel('Left/Right')
    pyplot.ylabel('In/Out')
    physics.set_axes_equal(ax)
    ax.set_zlabel('Up/Down')
    ax.plot(filtered_ball_position[0, :],
            filtered_ball_position[2, :], 'g--')

    pyplot.show()


def plot_all_vectors(time_vector,
                     position_vector,
                     velocity_vector):

    filtered_ball_position = position_vector[:, position_vector[1, :] > 0]
    filtered_ball_velocity = velocity_vector[:, position_vector[1, :] > 0]
    filtered_time = time_vector[position_vector[1, :] > 0]
    fig = pyplot.figure()
    ax = fig.add_subplot(121)

    ax.plot(filtered_time, filtered_ball_position[0, :], 'r', label='x')
    ax.plot(filtered_time, filtered_ball_position[1, :], 'b', label='y')
    ax.plot(filtered_time, filtered_ball_position[2, :], 'g', label='z')

    pyplot.grid()
    pyplot.xlabel('Time')
    pyplot.ylabel('Position (m)')
    pyplot.title('Position')
    pyplot.legend()

    ax = fig.add_subplot(122)
    ax.plot(filtered_time, filtered_ball_velocity[0, :], 'r', label='x')
    ax.plot(filtered_time, filtered_ball_velocity[1, :], 'b', label='y')
    ax.plot(filtered_time, filtered_ball_velocity[2, :], 'g', label='z')

    pyplot.grid()
    pyplot.xlabel('Time')
    pyplot.ylabel('Velocity (m/s)')
    pyplot.title('Velocity')

    pyplot.show()


class PositionObject:

    TIME_DELTA = .001

    def __init__(self,
                 position=numpy.zeros((3, 1)),
                 velocity=numpy.zeros((3, 1)),
                 acceleration=numpy.zeros((3, 1)),
                 initial_time=0,
                 time_delta=TIME_DELTA,
                 mass=10):

        self.position = position
        self.velocity = velocity
        self.acceleration = acceleration
        self.time = initial_time
        self.time_delta = time_delta
        self.mass = mass

    def propagate_position(self, delta_time):

        self.velocity += physics.calculate_gravity(self.mass) * delta_time
        self.position += self.velocity * delta_time


class Collector( PositionObject ):

    def __init__(self,
                 position=numpy.zeros((3, 1)),
                 velocity=numpy.zeros((3, 1)),
                 acceleration=numpy.zeros((3, 1)),
                 initial_time=0,
                 time_delta=PositionObject.TIME_DELTA ):

        super().__init__(position,
                         velocity,
                         acceleration,
                         initial_time,
                         time_delta)


class Target (PositionObject):

    def __init__(self,
                 position=numpy.zeros((3, 1)),
                 velocity=numpy.zeros((3, 1)),
                 acceleration=numpy.zeros((3, 1)),
                 initial_time=0,
                 time_delta=PositionObject.TIME_DELTA,
                 rcs_dBsm=0):

        super().__init__(position,
                         velocity,
                         acceleration,
                         initial_time,
                         time_delta)

        self.rcs_dBsm = rcs_dBsm


class GolfBall (Target):

    def __init__(self,
                 position=numpy.zeros((3, 1)),
                 velocity=numpy.zeros((3, 1)),
                 spin=numpy.zeros((3, 1)),
                 acceleration=numpy.zeros((3, 1)),
                 initial_time=0,
                 time_delta=Target.TIME_DELTA,
                 rcs_dBsm=0):

        super().__init__(position,
                         velocity,
                         acceleration,
                         initial_time,
                         time_delta,
                         rcs_dBsm)

        self.spin = spin

        self.length = .042  # meters
        self.surface_area = numpy.pi * (self.length / 2) ** 2
        self.ball_mass = .0433  # kg
        self.magnus_coefficient = .17
        self.drag_coefficient = .15  #

        # Needs to be moved to an environment class
        self.temperature = 21.  # celsius
        self.relative_humidity = .25  # fraction (0-1)
        self.pressure = 90097.16  # Pa

    def propagate_position(self, delta_time):

        air_density = physics.calculate_air_density(self.temperature,
                                                    self.relative_humidity,
                                                    self.pressure)

        air_resistance = -physics.calculate_air_resistance(self.drag_coefficient,
                                                           air_density,
                                                           self.surface_area,
                                                           self.velocity)

        magnus_force = physics.calculate_magnus_force(self.magnus_coefficient,
                                                      self.spin,
                                                      air_density,
                                                      self.surface_area,
                                                      self.velocity)

        resistance = air_resistance + magnus_force + physics.calculate_gravity(self.ball_mass)
        total_acceleration = resistance / self.ball_mass

        self.velocity[:, 0] += total_acceleration * delta_time
        self.position[:, 0] += self.velocity[:, 0] * delta_time

    def propogate_to_time(self, to_time):

        total_times = numpy.append(numpy.arange(self.time, to_time, self.time_delta), to_time)
        diff_times = numpy.diff(total_times)

        [self.propagate_position(diff_time) for diff_time in diff_times]

        print(self.position)






# Notes on ball flight:
#  Normal vector will be mostly perpendicular to initial velocity vector

def test_propogate_position():
    test=1

if __name__ == '__main__':

    print('Hello world')

    speed_of_light = 2.99875e8
    frequency_radiated = 2.0e9
    lambda_radiated = speed_of_light / frequency_radiated
    collector = Collector()
    launch_angle = 15 * numpy.pi / 180.
    launch_speed = 75
    target = GolfBall(numpy.array([[0.0, 0, 0]]).T,
                      numpy.array([[0.0, launch_speed*numpy.sin(launch_angle), launch_speed*numpy.cos(launch_angle)]]).T,
                      numpy.array([[-1, 0, 0]]).T)

    distance_between_positions = numpy.sqrt(numpy.sum(numpy.square(target.position) - numpy.square(collector.position)))
    phase_for_distance = 4 * numpy.pi * distance_between_positions / lambda_radiated

    sampling_period = .001

    time_to_analyze = numpy.arange(0, 10, sampling_period)
    time_deltas = numpy.diff(time_to_analyze)

    total_position = numpy.empty((3,len(time_to_analyze)))
    total_position[:, 0] = target.position[0, :]

    total_velocity = numpy.empty((3,len(time_to_analyze)))
    total_velocity[:, 0] = target.velocity[0, :]

    signal_return = numpy.zeros(time_to_analyze.shape, dtype=complex)

    for idx, diffs in enumerate(time_deltas):

        target.propagate_position(diffs)
        total_position[:, idx] = target.position[:, 0]
        total_velocity[:, idx] = target.velocity[:, 0]

        distance_between_positions = numpy.sqrt(numpy.sum(numpy.square(target.position) -
                                                          numpy.square(collector.position)))

        round_trip_time = 2 * distance_between_positions / speed_of_light
        phase_for_target = - 2 * numpy.pi * round_trip_time * frequency_radiated

        signal_return[idx] = 1 * numpy.exp(1j * phase_for_target)

        prf_hz = 1 / round_trip_time


    pyplot.figure()
    signal.spectrogram(signal_return)
    f, t, Sxx = signal.spectrogram(signal_return, 1. / sampling_period)
    compressedSignal = 20*numpy.log10(numpy.abs(Sxx))
    peakValue = numpy.max(compressedSignal, 0)
    peakIdx = numpy.argmax(compressedSignal, 0)
    f = numpy.fft.fftshift(f)
    frequency_detected = f[peakIdx]

    pyplot.pcolormesh(f, t, compressedSignal.T)
    pyplot.plot(frequency_detected, t,'r*')
    #pyplot.pcolormesh(t, f, compressedSignal)
    pyplot.xlabel('Frequency [Hz]')
    pyplot.ylabel('Time [sec]')

    pyplot.figure()
    pyplot.plot(numpy.linspace(-1/sampling_period/2, 1/sampling_period/2, 2048), 10*numpy.log10(numpy.abs(numpy.fft.fft(signal_return[0:2048]))))
    pyplot.grid()
    pyplot.axis('tight')

    plot_position_vector(total_position)
    plot_all_vectors(time_to_analyze,
                     total_position,
                     total_velocity)
    test=1
    pyplot.show()