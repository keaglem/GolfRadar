import numpy
import scipy
import matplotlib

# In this package, x is right (aft looking forward), y is up, and z is out

class PositionObject:

    TIME_DELTA = .001

    def __init__(self,
                 position=numpy.zeros((3, 1)),
                 velocity=numpy.zeros((3, 1)),
                 acceleration=numpy.zeros((3, 1)),
                 initial_time=0,
                 time_delta=TIME_DELTA):

        self.position = position
        self.velocity = velocity
        self.acceleration = acceleration
        self.time = initial_time
        self.time_delta = time_delta

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


speed_of_light = 2.99875e8
frequency_radiated = 16.0e9
lambda_radiated = speed_of_light / frequency_radiated
collector = Collector()
target = Target(numpy.array([[0, 20, 50]]).T,
                numpy.array([[0, 1, 5]]).T)

distance_between_positions = numpy.sqrt(numpy.sum(numpy.square(target.position) - numpy.square(collector.position)))
phase_for_distance = 4 * numpy.pi * distance_between_positions / lambda_radiated



# Notes on ball flight:
#  Normal vector will be mostly perpendicular to initial velocity vector

def test_propogate_position():
    test=1

if __name__ == '__main__':

    print('Hello world')


