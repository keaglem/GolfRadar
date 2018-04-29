import numpy
import scipy
import matplotlib
from physics import physics

# In this package, x is right (aft looking forward), y is up, and z is out

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
                 acceleration=numpy.zeros((3, 1)),
                 initial_time=0,
                 time_delta=Target.TIME_DELTA,
                 rcs_dBsm=0,
                 spin=numpy.zeros((3,1))):

        super().__init__(position,
                         velocity,
                         acceleration,
                         initial_time,
                         time_delta,
                         rcs_dBsm)

        self.spin = spin

        self.length = .042  # meters
        self.surface_area = 4 * numpy.pi * (self.length / 2) ** 2
        self.ball_mass = .0433  # kg
        self.magnus_coefficient = .17
        self.drag_coefficient = .2  #

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

        resistance = air_resistance + magnus_force + physics.calculate_gravity(self.mass)
        total_acceleration = resistance / self.mass

        self.velocity[:, 0] += total_acceleration * delta_time
        self.position[:, 0] += self.velocity[:, 0] * delta_time

    def propogate_to_time(self, to_time):

        total_times = numpy.append(numpy.arange(self.time, to_time, self.time_delta), to_time)
        diff_times = numpy.diff(total_times)

        [self.propagate_position(diff_time) for diff_time in diff_times]

        print(self.position)


speed_of_light = 2.99875e8
frequency_radiated = 16.0e9
lambda_radiated = speed_of_light / frequency_radiated
collector = Collector()
target = GolfBall(numpy.array([[0.0, 20, 50]]).T,
                  numpy.array([[0.0, 0, 50]]).T)

distance_between_positions = numpy.sqrt(numpy.sum(numpy.square(target.position) - numpy.square(collector.position)))
phase_for_distance = 4 * numpy.pi * distance_between_positions / lambda_radiated

target.propogate_to_time(2)

test=1


# Notes on ball flight:
#  Normal vector will be mostly perpendicular to initial velocity vector

def test_propogate_position():
    test=1

if __name__ == '__main__':

    print('Hello world')


