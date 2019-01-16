class Component(object):
    '''Container for forces in an aircraft component'''

    def __init__(self, name):
        '''Constructor'''
        self.box = None
        self.zonelist = []
        self.name = name
        self.clear()

    def clear(self):
        '''Clear forces contained in the component'''
        self.force_friction = [0, 0, 0]
        self.force_pressure = [0, 0, 0]
        self.area = 0
        self.moment = 0.

    def total_force(self, vec):
        '''Compute total force in component along vector

        arguments:
            vec (list of float): vector (lift or drag)
        '''
        return self.friction_force(vec) + self.pressure_force(vec)

    def friction_force(self, vec):
        '''Compute friction force along vector'''
        return (self.force_friction[0]*vec[0] +
                self.force_friction[1]*vec[1] +
                self.force_friction[2]*vec[2]
                )

    def pressure_force(self, vec):
        '''Compute pressure force along vector'''
        return (self.force_pressure[0]*vec[0] +
                self.force_pressure[1]*vec[1] +
                self.force_pressure[2]*vec[2]
                )
