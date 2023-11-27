# Default solar system objects (Obviously incomplete)

from simulator import AstronomicalObject
from scipy import constants
import numpy as np

sun = AstronomicalObject("Sun",
                            mass=1.989e+30,
                            position=np.array([0., 0]),
                            velocity=np.array([0., 0]),
                            radius=696_340_000,
                            colour = 'y',
                            colour2 = 'orange'
                        )

earth = AstronomicalObject("Earth",
                           mass=5.972e+24,
                           position=np.array([constants.au, 0.]),
                           velocity=np.array([0., 107_226*constants.kmh]),
                           radius=None,
                           colour = 'b',
                           colour2 = 'green'
                        )

moon = AstronomicalObject("Moon",
                            mass=7.34767309e+22,
                            position=np.array([constants.au + 384_400_000, 0.]),
                            velocity=np.array([0., 107_226*constants.kmh + 3683*constants.kmh]),
                            radius=None,
                            colour = 'grey',
                            colour2 = 'grey'
                        )
