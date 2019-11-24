import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants:
G = 6.67e-11
pc = 3.086e16  # Parsec
sm = 1.989e30  # Solar Mass
yr = 60 * 60 * 24 * 365  # 1 year in seconds

# Simulation time conditions:
time_step = 2e5 * yr  # Time between each step.
time_run = 2e9 * yr  # Total time simulation is run for.
frames = 10  # Number of intervals between images being shown.
no_step = time_run / time_step  # Total number of steps in simulation.
interval = no_step / frames  #
image_time_step = time_run / (frames * 1e9 * yr)  # Time between images being shown.

# Galaxy 1 starting conditions:
mg1 = 1e11 * sm  # Mass of primary galaxy.
xg1 = 0  # x position of primary galaxy.
yg1 = 0  # y position of primary galaxy.
zg1 = 0  # y position of primary galaxy.
vxg1 = 0  # x velcoity of primary galaxy.
vyg1 = 75e3  # y velcotiy of primary galaxy.
vzg1 = 0  # z velocity of primary galaxy.

# Galaxy 2 starting conditions:
mg2 = 1e11 * sm  # Mass of secondary galaxy.
xg2 = 60e3 * pc  # x position of secondary galaxy.
yg2 = 250e3 * pc  # y position of secondary galaxy.
zg2 = 0  # z position of secondary galaxy.
vxg2 = 0  # x velcoity of secondary galaxy.
vyg2 = -75e3  # y velocity of secondary galaxy.
vzg2 = 0  # z velocity of secondary galaxy.

# Galaxy ring conditions
no_rings = 5  # Number of rings.
no_rp = 6  # Number of particles in innermost ring.
ring_rad = 2.5e3 * pc  # Radius of innermost ring from galaxy centre.

objects = []  # List of all objects in galaxy.


class Body:
    def __init__(self, name, m, x, y, z, vx, vy, vz, saved_x, saved_y, saved_z, image_x, image_y, image_z, colour):
        self.name = name  # Name of body.
        self.m = m  # Mass of body.
        self.x = x  # x position of body.
        self.y = y  # y position of body.
        self.z = z  # z position of body.
        self.vx = vx  # x velocity of body.
        self.vy = vy  # y velocity of body.
        self.vz = vz  # z velocity of body.
        self.saved_x = []  # Array of all x position values of body.
        self.saved_y = []  # Array of all y position values of body.
        self.saved_z = []  # Array of all z position values of body.
        self.image_x = []  # Array of all x position values of body, used for images.
        self.image_y = []  # Array of all y position values of body, used for images.
        self.image_z = []  # Array of all z position values of body, used for images.
        self.colour = colour  # Colour of body on images.

    def name(self):
        return self.name

    def mass(self):
        return self.m

    def x(self):
        return self.x

    def y(self):
        return self.y

    def z(self):
        return self.z

    def vx(self):
        return self.vx

    def vy(self):
        return self.vy

    def vz(self):
        return self.vz

    def saved_x(self):
        return self.saved_x

    def saved_y(self):
        return self.saved_y

    def saved_z(self):
        return self.saved_z

    def image_x(self):
        return self.image_x

    def image_y(self):
        return self.image_y

    def image_z(self):
        return self.image_z

    def colour(self):
        return self.colour

    def force_calc(self, other):  # Force calculation between two bodies
        rx = self.x - other.x
        ry = self.y - other.y  # Distance between two bodies in all directions.
        rz = self.z - other.z
        f = -(G * self.m * other.m)/(rx ** 2 + ry ** 2 + rz ** 2)  # Total force calculation.
        theta = math.atan2(ry, rx)  # Azimuthal angle.
        phi = math.acos(rz / math.sqrt(rx ** 2 + ry ** 2 + rz ** 2))  # Polar angle.
        fx = math.cos(theta) * math.sin(phi) * f
        fy = math.sin(theta) * math.sin(phi) * f  # Force calculation in each direction.
        fz = math.cos(phi) * f
        return fx, fy, fz


# class Cell:
#     def __init__(self, min, width, pp, pm):
#         self.min = []
#         self.width = width
#         self.centre = self.centre_calc
#         self.no_part = 0
#         self.pp = [[]]
#         self.pm = []
#         self.ppm = self.ppm_calc
#         self.M = np.sum(self.pm)
#         self.com = self.com_calc
#         self.children = []
#         self.leaf = False
#
#     def min(self):
#         return self.min
#
#     def width(self):
#         return self.width
#
#     def centre(self):
#         return self.centre
#
#     def no_part(self):
#         return self.no_part
#
#     def pm(self):
#         return self.pm
#
#     def M(self):
#         return self.M
#
#     def com(self):
#         return self.com
#
#     def children(self):
#         return self.children
#
#     def leaf(self):
#         return self.leaf
#
#     def centre_calc(self):
#         self.centre[0] = self.min[0] + (self.width / 2)
#         self.centre[1] = self.min[1] + (self.width / 2)
#         self.centre[2] = self.min[2] + (self.width / 2)
#
#     def ppm_calc(self):
#         for j in range(len(pm)):
#             self.ppm[0] = pp[j][0] * pm[j]
#             self.ppm[1] = pp[j][1] * pm[j]
#             self.ppm[2] = pp[j][2] * pm[j]
#
#     def com_calc(self):
#         self.com[0] = self.ppm[0] / self.M
#         self.com[1] = self.ppm[1] / self.M
#         self.com[2] = self.ppm[2] / self.M
#
#     def create_children(self):
#         aaaaaaaaaaaaaaaaaaaaaa


def create_galaxies():  # Creates two interacting galaxies using Body class.
    objects.append(Body("Primary", mg1, xg1, yg1, zg1, vxg1, vyg1, vzg1, 0, 0, 0, 0, 0, 0, 'bo'))
    # Creates primary galaxy and adds it to list of bodies.
    objects.append(Body("Secondary", mg2, xg2, yg2, zg2, vxg2, vyg2, vzg2, 0, 0, 0, 0, 0, 0, 'ro'))
    # Creates secondary galaxy and adds it to list of bodies.


def create_rings():  # Creates set of rings for each galaxy.
    for k in range(0, no_rings):  # Creating rings for primary galaxy.
        r = (k + 1) * ring_rad  # Radius of each ring.
        n = (k + 1) * no_rp  # Number of particles in each ring.
        v = math.sqrt(G * mg1 / r)  # Velocity of particles in each ring.
        for i in range(0, n):
            angle = 2 * math.pi * i / n  # Assigning particles to a ring formation.
            x = r * math.cos(angle)
            y = r * math.sin(angle)  # x and y coordinates of each particle in ring.
            vx = v * math.sin(angle)
            vy = -v * math.cos(angle)  # x and y velocities of each particle in ring.
            objects.append(Body("Test", 1, xg1 + x, yg1 + y, zg1, vxg1 + vx, vyg1 + vy, vzg1, 0, 0, 0, 0, 0, 0, 'g.'))
                                                     # Adding particles to list of Bodies. Test particles, mass = 1kg.

    # for k in range(0, no_rings):  # Creating rings for secondary galaxy.
    #     r = (k + 1) * ring_rad  # Radius of each ring.
    #     n = (k + 1) * no_rp  # Number of particles in each ring.
    #     v = math.sqrt(G * mg2 / r)  # Velocity of particles in each ring.
    #     for i in range(0, n):
    #         angle = 2 * math.pi * i / n  # Assigning particles to a ring formation.
    #         x = r * math.cos(angle)
    #         y = r * math.sin(angle)  # x and y coordinates of each particle in ring.
    #         vx = v * math.sin(angle)
    #         vy = -v * math.cos(angle)  # x and y velocities of each particle in ring.
    #         objects.append(Body('Test', 1, xg2 + x, yg2 + y, zg2, vxg2 + vx, vyg2 + vy, vzg2, 0, 0, 0, 0, 0, 0, 'yo'))
                                                    # Adding particles to list of Bodies. Test particles, mass = 1kg.


def leapfrog_initial(bodies):  # Produces kick start for the leapfrog algorithm.
    print("Calculating...")
    for body in bodies:
        body.saved_x.append(body.x)
        body.saved_y.append(body.y)  # Appends initial positions of each body to a
        body.saved_z.append(body.z)  # list of saved positions for that body.
        body.image_x.append(body.x)
        body.image_y.append(body.y)  # Appends initial positions of each body to a
        body.image_z.append(body.z)  # list of plotted positions for that body.

    force = {}
    for body in bodies:
        total_fx = total_fy = total_fz = 0.0  # Sets force for all bodies to 0.
        for other in bodies:
            if body is other:  # Checking that not calculating force of a body on itself.
                continue
            if other.name != "Test":  # Does not calculate force due to ring/test particles.
                fx, fy, fz = body.force_calc(other)
                total_fx += fx
                total_fy += fy  # Add together forces of al other bodies acting on that body.
                total_fz += fz
        force[body] = (total_fx, total_fy, total_fz)

    for body in bodies:  # Kick start position and velocities for all bodies.
        fx, fy, fz = force[body]

        body.vx += (fx / body.m) * (time_step / 2)
        body.vy += (fy / body.m) * (time_step / 2)  # Calculates initial half step for velocity in each direction,
        body.vz += (fz / body.m) * (time_step / 2)  # using v_0.5 = v_0 + a * (dt/2), where a = F(x_0)/m.

        body.x += body.vx * time_step
        body.y += body.vy * time_step  # Uses half step in velocity to calculate new position,
        body.z += body.vz * time_step  # using x_1 = x_0 + v_0.5 * dt.

        body.x = body.x - objects[0].saved_x[-1]
        body.y = body.y - objects[0].saved_y[-1]
        body.z = body.z - objects[0].saved_z[-1]

        body.saved_x.append(body.x)
        body.saved_y.append(body.y)  # Saving new position to list of previous positions of body.
        body.saved_z.append(body.z)


def leapfrog(bodies):  # Updates the position and velocity of each particle using a leapfrog algorithm.
    step = 1
    percent = 0.0
    print(percent)
    while True:
        step += 1
        for a in range(100):
            if step == (a * no_step) / 100:
                percent += 1.0
                print(percent)  # Percentage calculator for running the simulation.
        if step == no_step:
            return  # Stop simulation when all steps done.
        else:
            force = {}
            for body in bodies:
                total_fx = total_fy = total_fz = 0.0  # Sets force for all bodies to 0.
                for other in bodies:
                    if body is other:  # Checking that not calculating force of a body on itself.
                        continue
                    if other.name != "Test":  # Does not calculate force due to ring/test particles.
                        fx, fy, fz = body.force_calc(other)
                        total_fx += fx
                        total_fy += fy  # Add together forces of al other bodies acting on that body.
                        total_fz += fz
                force[body] = (total_fx, total_fy, total_fz)

            for body in bodies:
                fx, fy, fz = force[body]

                body.vx += (fx / body.m) * time_step
                body.vy += (fy / body.m) * time_step  # Calculates the new velocity in each direction,
                body.vz += (fz / body.m) * time_step  # using v_i+1.5 = v_i+0.5 + a * dt, where a = F(x_i)/m.

                body.x += body.vx * time_step
                body.y += body.vy * time_step  # Uses new velocity to calculate new position,
                body.z += body.vz * time_step  # using x_i+1 = x_i + v_i+0.5 * dt.

                body.x = body.x - objects[0].saved_x[-1]
                body.y = body.y - objects[0].saved_y[-1]  # THIS WORKS THE SAME AS IN NBCode3.py !!!!!!! SIMPLER!!!
                body.z = body.z - objects[0].saved_z[-1]

                body.saved_x.append(body.x)
                body.saved_y.append(body.y)  # Saving new position to list of previous positions of body.
                body.saved_z.append(body.z)

                if step % interval == 0 or step == no_step - 1:
                    body.image_x.append(body.x)
                    body.image_y.append(body.y)  # Saving new position to list of positions that will be plotted.
                    body.image_z.append(body.z)


def pericentre_calc():  # Calculates the pericentre of the interaction and the time at which it occurs.
    pericentres = []  # Create list of pericentres at each time step.
    for i in range(len(objects[0].saved_x)):
        xpc = (objects[0].saved_x[i] - objects[1].saved_x[i]) ** 2
        ypc = (objects[0].saved_y[i] - objects[1].saved_y[i]) ** 2
        zpc = (objects[0].saved_z[i] - objects[1].saved_z[i]) ** 2
        pericentre = math.sqrt(xpc + ypc + zpc)  # Calculate pericentre.
        pericentres.append(pericentre)
    min_pericentre = min(pericentres)  # Minimum value of pericentre found in list.
    min_pericentre_pc = min_pericentre / (1e3 * pc)
    t = ((pericentres.index(min(pericentres)) + 1) * time_step) / (1e9 * yr)  # Min pericentre position in list
                                                                              # used to find time.

    print("\nPericentre =", round(min_pericentre_pc, 2), "kpc at t =", round(t, 2), "Gyrs.\n")

    print("The time between images is", image_time_step, "Gyrs\n")  # Prints time in between images.


def plot():  # Plot images of interaction.
    print("Producing Images.")
    for j in range(len(objects)):  # Shows last image first, so plot in reverse to see in order.
        objects[j].saved_x.reverse()
        objects[j].saved_y.reverse()
        objects[j].saved_z.reverse()
        objects[j].image_x.reverse()
        objects[j].image_y.reverse()
        objects[j].image_z.reverse()
    for i in range(len(objects[0].image_x)):  # Plots data of each time step on a separate figure.
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        title = ((-i + 10) * image_time_step)
        ax.set_title('t = %.2f Gyrs' % title, fontsize=10)  # Title the figure with the time of interaction.
        ax.set_xlabel('X', fontsize=10)
        ax.set_ylabel('Y', fontsize=10)  # Set axis labels.
        ax.set_zlabel('Z', fontsize=10)
        ax.set_xlim(-2e21, 2e21)
        ax.set_ylim(-2e21, 2e21)  # Set axis limits.
        ax.set_zlim(-2e11, 2e11)
        ax.plot3D(objects[0].saved_x, objects[0].saved_y, objects[0].saved_z, 'b-')  # Plotting path of primary galaxy.
        ax.plot3D(objects[1].saved_x, objects[1].saved_y, objects[1].saved_z, 'r-')  # Plotting path of secondary
                                                                                     # galaxy.
        for j in range(len(objects)):  # Plotting all objects on figure.
            a = [objects[j].image_x[i]]
            b = [objects[j].image_y[i]]
            c = [objects[j].image_z[i]]
            ax.plot3D(a, b, c, objects[j].colour)


def main():  # Calling all functions in order.

    create_galaxies()

    create_rings()

    leapfrog_initial(objects)

    leapfrog(objects)

    pericentre_calc()

    plot()

    plt.show()


if __name__ == '__main__':
    main()