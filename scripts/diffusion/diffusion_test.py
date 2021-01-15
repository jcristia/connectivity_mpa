# compare horizontal_diffusion and current_uncertainty and decide which one to use

# horizontal diffusion was added September 2020:
# https://github.com/OpenDrift/opendrift/commit/10fed704d2c45f1505b365752d53ed683220021e#diff-23e8c0e5a71f7b65333678fdd6884b5270d3714e294ce1b5c6eaaeb2f9ef3462
# " * New config setting `drift:horizontal_diffusivity`, providing time-step independent diffusion, in contrast to `drift:current_uncertainty` and `drift:wind_uncertainty` "
# I don't understand why it says "time-step independent" since it also uses the time step in its calculation.

# Aside from the differences I detail below, they are different in the where they get called.
# they both get called within run(). Within run() there is a for loop for each time step.

### they both end up eventually calling update_position() in basemodel
# horizontal_diffusion:
# horizontal_diffusion gets called in the for loop for each time step
# it calculates x_vel and y_vel for JUST movement by diffusion and sends it to update_positions()
# current uncertainty:
# in the for loop, get_environment gets called, which adds the random kick value to env[x_sea_water_velocity], so that it is actually part of the variable now, which it returns
# then back in the for loop, update() gets called, which is in oceandrift.py
# update calls advect_ocean_current(), which is in physics_methods.py
# advect_ocean_current() does the runge kutta method. To calculate midpoints it
# calls get_environment again, but it sets profiles=None so that current_uncertainty isn't calculated each time
# it then calls update_positions back in basemodel with the rk averaged midpoint velocities

# SO THE KEY DIFFERENCE IS (aside from the way they pull from the random distribution):
# CU first adds the random kick in before RK is calculated so that the overall speed and direction is changed before the new position is calculated.
# HD adds it in after RK is done. So we first find where it moves based on just advection, the it adds in a small random kick.

# I feel like current_uncertainty is aptly named because it changes u and v
# Therefore, HD might be the better way to go simply because it is more intuitive.
# However, I'm sure either way is legit.



# code from Opendrift basemodel for horizontal diffusion
# the equation is exactly the same as what Ben provided me
# However, it takes the values in m/2 and multiplies it by a random number with a standard deviation of 1
# it the calls update_positions
def horizontal_diffusion(self):
    """Move elements with random walk according to given horizontal diffuivity."""
    D = self.get_config('drift:horizontal_diffusivity')
    if D == 0:
        self.logger.debug('Horizontal diffusivity is 0, no random walk.')
        return
    dt = self.time_step.total_seconds()
    x_vel = np.sqrt(2*D/dt)*np.random.normal(scale=1, size=self.num_elements_active())
    y_vel = np.sqrt(2*D/dt)*np.random.normal(scale=1, size=self.num_elements_active())
    speed = np.sqrt(x_vel*x_vel+y_vel*y_vel)
    self.logger.debug('Moving elements according to horizontal diffusivity of %s, with speeds between %s and %s m/s'
                        % (D, speed.min(), speed.max()))
    self.update_positions(x_vel, y_vel)


# code from Opendrift basemodel for current uncertainty
# it takes the provided value in m/s, uses it as the standard deviation in a distribution
if 'x_sea_water_velocity' in variables and \
        'y_sea_water_velocity' in variables:
    std = self.get_config('drift:current_uncertainty')
    if std > 0:
        self.logger.debug('Adding uncertainty for current: %s m/s' % std)
        env['x_sea_water_velocity'] += np.random.normal(
            0, std, self.num_elements_active())
        env['y_sea_water_velocity'] += np.random.normal(
            0, std, self.num_elements_active())
    std = self.get_config('drift:current_uncertainty_uniform')
    if std > 0:
        self.logger.debug('Adding uncertainty for current: %s m/s' % std)
        env['x_sea_water_velocity'] += np.random.uniform(
            -std, std, self.num_elements_active())
        env['y_sea_water_velocity'] += np.random.uniform(
            -std, std, self.num_elements_active())


# use the two methods above to calculate two distributions and compare
# What I need to investigate:
#   horiz diff does 0.25 * a random number with a std of 1. Whereas current_uncertainty does just a random number with a std of 0.25. How different are these?
# what are the max values for each distrubition? 3 standard devs is 99% of the mean, so for horiz diff, that is 3, for CU that is 0.75.
# Perhaps even though you end up selecting MORE big numbers with a std of 1, once they are reduced by multiplying by 0.25, the overall distributions comes out the same.
# SO THE FIRST TEST SHOULD JUST BE selecting a certain amount with both methods and comparing their distributions to see if they are similar.

import numpy as np
D = 1.5
dt = 60
vel = np.sqrt(2*D/dt)
HD = np.sqrt(2*D/dt)*np.random.normal(scale=1, size=1000)
CU = np.random.normal(0, np.sqrt(2*D/dt), size=1000)

import matplotlib.pyplot as plt
plt.hist(HD, bins=20)
plt.show
plt.hist(CU, bins=20)
plt.show

# run the above code 3-4 times
# It looks like these aren't that different

# DECISION: I will go with horizontal_diffusion because:
# I like the way it gets calculated in opendrift at the very end. This is more what I think of when I think of a random kick.
# Also, it will be easier to write about than explaining how I calculated current_uncertainty.

# I compared the two in a simulation and they are pretty much identical.