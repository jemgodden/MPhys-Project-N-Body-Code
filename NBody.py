import math
import numpy as np
import matplotlib.pyplot as plt
import time
import os
from mpl_toolkits.mplot3d import Axes3D
from IntData import *

objects = []  # List of all objects in galaxy.

start_time = time.time()  # Sets start time in order to find runtime of program.


class Body:
    def __init__(self, name, m, position, velocity, colour):
        self.name = name  # Name of body.
        self.m = m  # Mass of body.
        self.xyz = position  # Array of x, y and z position of body.
        self.vxyz = velocity  # Array of x, y and z velocities of body.
        self.saved_xyz = [[], [], []]  # Array of all x, y and z position values of body.
        self.force = [0.0] * 3  # A list of all the forces, in each direction, for the body at each time step.
        self.colour = colour  # Colour of body on images.

    def name(self):
        return self.name

    def mass(self):
        return self.m

    def xyz(self):
        return self.xyz

    def vxyz(self):
        return self.vxyz

    def saved_xyz(self):
        return self.saved_xyz

    def force(self):
        return self.force

    def colour(self):
        return self.colour

    def force_calc(self, other):  # Force calculation between two bodies
        fxyz = [None] * 3
        rx = self.xyz[0] - other.xyz[0]
        ry = self.xyz[1] - other.xyz[1]  # Distance between two bodies in all directions.
        rz = self.xyz[2] - other.xyz[2]
        r = rx ** 2 + ry ** 2 + rz ** 2
        f = -(G * self.m * other.m) / r  # Total force calculation.
        theta = math.atan2(ry, rx)  # Azimuthal angle.
        phi = math.acos(rz / math.sqrt(r))  # Polar angle.
        fxyz[0] = math.cos(theta) * math.sin(phi) * f
        fxyz[1] = math.sin(theta) * math.sin(phi) * f  # Force calculation in each direction.
        fxyz[2] = math.cos(phi) * f

        for i in range(len(fxyz)):  # Checks whether each particle is a test and adds forces where needed.
            if self.name != "Test" and other.name != "Test":
                self.force[i] += fxyz[i]
                other.force[i] -= fxyz[i]
            if self.name != "Test" and other.name == "Test":
                other.force[i] -= fxyz[i]
            if self.name == "Test" and other.name != "Test":
                self.force[i] += fxyz[i]


def make_directories():  # Checks directories exist and makes them if not.
    if not os.path.exists('./Forwards'):
        os.makedirs('./Forwards')
    if not os.path.exists('./Backwards'):
        os.makedirs('./Backwards')


def create_galaxies():  # Creates two interacting galaxies using Body class.
    objects.append(Body("Primary", mg1, [xg1, yg1, zg1], [vxg1, vyg1, vzg1], 'bo'))
    # Creates primary galaxy and adds it to list of bodies.
    objects.append(Body("Secondary", mg2, [xg2, yg2, zg2], [vxg2, vyg2, vzg2], 'ro'))
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
            objects.append(Body("Test", 1, [xg1 + x, yg1 + y, zg1], [vxg1 + vx, vyg1 + vy, vzg1], 'c.'))
                                                # Adding ring particles to list of Bodies. Test particles, mass = 1kg.

    if secondary_disk:
        for k in range(0, no_rings):  # Creating rings for secondary galaxy.
            r = (k + 1) * ring_rad  # Radius of each ring.
            n = (k + 1) * no_rp  # Number of particles in each ring.
            v = math.sqrt(G * mg2 / r)  # Velocity of particles in each ring.
            for i in range(0, n):
                angle = 2 * math.pi * i / n  # Assigning particles to a ring formation.
                x = r * math.cos(angle)
                y = r * math.sin(angle)  # x and y coordinates of each particle in ring.
                vx = v * math.sin(angle)
                vy = -v * math.cos(angle)  # x and y velocities of each particle in ring.
                objects.append(Body('Test', 1, [xg2 + x, yg2 + y, zg2], [vxg2 + vx, vyg2 + vy, vzg2], 'm.'))
                                                # Adding ring particles to list of Bodies. Test particles, mass = 1kg.


def rewind_initial():  # Automatically creates the file with initial conditions to run a simulation backwards.
    open("Initial_Conditions.txt", "w+").close()  # Clears contents of file.

    file_name = - time_run / Gyr
    file1 = open("Forwards/image_%.2f.txt" % file_name)  # Reads file at end of last forward simulation.
    lines = file1.readlines()
    file2 = open("Initial_Conditions.txt", "w")
    file2.writelines(lines)  # Writes file full of information for backwards interaction.
    file1.close()
    file2.close()

    read_initial_conditions()


def read_initial_conditions():  # Reads in initial conditions, of all particles, to the simulation form a text file.
    file = open("Initial_Conditions.txt", "r")
    for line in file:  # Reads file line by line.
        data = line.strip().split()
        objects.append(Body(data[0], float(data[1]), [float(data[2]), float(data[3]), float(data[4])],
                            [float(data[5]), float(data[6]), float(data[7])], data[8]))  # Creates bodies using
    file.close()                                                                         # information read from file.


def leapfrog_initial(bodies, step):  # Produces kick start for the leapfrog algorithm.
    print("Calculating...")
    position_print(bodies, step)

    for body in bodies:
        body.saved_xyz[0].append(body.xyz[0])
        body.saved_xyz[1].append(body.xyz[1])  # Appends initial positions of each body to a
        body.saved_xyz[2].append(body.xyz[2])  # list of saved positions for that body.

    for i in range(len(objects) - 1):  # Loop that does not calculate each force more than once.
        body = objects[i]
        for j in range(i + 1, len(objects)):
            other = objects[j]
            if body.name == "Test" and other.name == "Test":  # Does not calculate force is both particles are test.
                continue
            else:
                body.force_calc(other)

    for body in bodies:  # Kick start position and velocities for all bodies.
        for i in range(len(body.force)):
            body.vxyz[i] += (body.force[i] / body.m) * (time_step / 2)  # Calculate initial half-step velocity in each direction.
            body.xyz[i] += body.vxyz[i] * time_step  # Uses half step in velocity to calculate new position.
            body.saved_xyz[i].append(body.xyz[i])  # Saving new position to list of previous positions of body.
        for j in range(len(body.force)):
            body.force[j] = 0


def leapfrog(bodies):  # Updates the position and velocity of each particle using a leapfrog algorithm.
    step = 0
    leapfrog_initial(bodies, step)

    percent = 0.0
    print(percent)
    while True:
        step += 1
        for a in range(100):
            if step == (a * no_step) / 100:
                percent += 1.0
                print(percent)  # Percentage calculator for running the simulation.

        if step == no_step:
            position_print(bodies, step)
            return  # Stop simulation when all steps done.
        else:
            for i in range(len(objects) - 1):  # Loop that does not calculate each force more than once.
                body = objects[i]
                for j in range(i + 1, len(objects)):
                    other = objects[j]
                    if body.name == "Test" and other.name == "Test":  # Does not calculate force of two test particles.
                        continue
                    else:
                        body.force_calc(other)

            for body in bodies:
                for i in range(len(body.force)):
                    body.vxyz[i] += (body.force[i] / body.m) * time_step  # Calculates new velocity in each direction.
                    body.xyz[i] += body.vxyz[i] * time_step  # Uses new velocity to calculate new position.
                    body.saved_xyz[i].append(body.xyz[i])  # Saving new position to list of previous positions of body.
                for j in range(len(body.force)):
                    body.force[j] = 0

            if step % int(interval) == 0:
                position_print(bodies, step)  # Print information on particles to a file at a particular time.


def info():
    x = pericentre_calc()
    min_pericentre_pc = x[0]
    t = x[1]

    total_time = (time.time() - start_time) / 60  # Calculates runtime in minutes.

    print("\nThe runtime for this simulation was %.2f minutes.\n\n" % total_time)  # Prints runtime.

    print("There are", tot_rp, "particles in one galaxy's disk.\n")  # Prints total number of particles in a
                                                                         # galaxy's disk.
    print("There are a total of", tot_part, "particles in the simulation.\n")  # Prints total number of particles
                                                                               # in a simulation.
    print("Pericentre of interaction =", round(min_pericentre_pc, 2), "kpc, occurring at t =", round(t, 2), "Gyrs.\n")
                                                            # Prints value of pericentre and time at which it occurs.
    print("The time between each of the", (frames + 1), "images is", image_time_step, "Gyrs.\n\n")
                                                            # Prints number of images and time between them.
    if rewind:
        file = open("Backwards/RewindPericentreInfo.txt", "w+")  # Writes information to file for backwards interaction.
    else:
        file = open("Forwards/PericentreInfo.txt", "w+")  # Writes information to file for forwards interaction.
    file.write("{0} {1} {2} {3}".format(tot_rp, tot_part, min_pericentre_pc, t))  # Prints information values to a file.
    file.close()


def pericentre_calc():  # Calculates the pericentre of the interaction and the time at which it occurs.
    pericentres = []  # Create list of pericentres at each time step.
    xyz_pc = [None] * 3
    for i in range(len(objects[0].saved_xyz[0])):
        for j in range(len(xyz_pc)):
            xyz_pc[j] = (objects[0].saved_xyz[j][i] - objects[1].saved_xyz[j][i]) ** 2
        pericentre = math.sqrt(xyz_pc[0] + xyz_pc[1] + xyz_pc[2])  # Calculates pericentre.
        pericentres.append(pericentre)

    min_pericentre = min(pericentres)  # Minimum value of pericentre found in list.
    min_pericentre_pc = min_pericentre / kpc
    t = ((pericentres.index(min(pericentres)) + 1) * time_step) / Gyr  # Min pericentre position in list used
                                                                       # to find time.
    return min_pericentre_pc, t


def path_print():
    if rewind:
        file1 = open("Backwards/RWPriGalPath.txt", "w+")
    else:
        file1 = open("Forwards/PriGalPath.txt",
                     "w+")  # Prints coordinates of primary galaxy at every time step to a file.
    for i in range(len(objects[0].saved_xyz[0])):
        file1.write("{0} {1} {2}\n".format(objects[0].saved_xyz[0][i], objects[0].saved_xyz[1][i],
                                           objects[0].saved_xyz[2][i]))
    file1.close()

    if rewind:
        file2 = open("Backwards/RWSecGalPath.txt", "w+")
    else:
        file2 = open("Forwards/SecGalPath.txt",
                     "w+")  # Prints coordinates of secondary galaxy at every time step to a file.
    for i in range(len(objects[1].saved_xyz[0])):
        file2.write("{0} {1} {2}\n".format(objects[1].saved_xyz[0][i], objects[1].saved_xyz[1][i],
                                           objects[1].saved_xyz[2][i]))
    file2.close()


def position_print(bodies, step):  # Printing information of each particle at the time an image is seen.
    txt_title = step * time_step / Gyr
    if rewind:
        file = open("Backwards/rimage_%.2f.txt" % txt_title, "w+")
    else:
        file = open("Forwards/image_%.2f.txt" % txt_title, "w+")  # Opens/creates a file with the name of the time of image.
    for body in bodies:
        file.write("{0} {1} {2} {3} {4} {5} {6} {7} {8}\n".format(body.name, body.m, body.xyz[0], body.xyz[1],
                                            body.xyz[2], body.vxyz[0], body.vxyz[1], body.vxyz[2], body.colour))
    file.close()                    # Writes all information on a particle to one line.


def main():  # Calling all functions in order.

    make_directories()

    if initial_txt:
        read_initial_conditions()
    if rewind:
        rewind_initial()
    else:
        create_galaxies()

        create_rings()

    leapfrog(objects)

    info()

    path_print()

    print("The data has been printed to text files. Please use the Plotter function to see the images.\n")


if __name__ == '__main__':
    main()
