# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib as mpl
from matplotlib.backends.backend_pgf import FigureCanvasPgf
mpl.backend_bases.register_backend('pdf', FigureCanvasPgf)
import matplotlib.pyplot as plt
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['text.usetex'] = True
mpl.rcParams['pgf.texsystem'] = 'lualatex'
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams.update({'figure.autolayout': True, 'legend.fontsize': 20})


class data_processing():
    def __init__(self, root, tau):
        self.root = root
        self.tau = tau
        self.time = []
        self.time_nd = []
        self.mass = []
        self.num_particles = 0

    def sort_files(self):
        def key_func(x):
            return int(x.strip('.txt').split('_')[-1])

        self.files = [f for f in os.listdir(self.root)
                      if os.path.isfile(os.path.join(self.root, f))]
        if any("setup" in filename for filename in self.files):
            a = [x for x in self.files if "setup" in x]
        try:
            self.files.remove(a[0])
        except ValueError:
            print("No setup file")

        self.files.sort(key=key_func)  # sort files into "natural order"

    def read_files(self):
        for file in self.files:
            if file.endswith('.txt'):
                with open(os.path.join(self.root, file), 'r') as f:
                    self.time.append(int(file.strip('.txt')
                                     .split('_')[-1])/1000000)
                    self.time_nd.append(int(file.strip('.txt')
                                        .split('_')[-1])/1000000/self.tau)
                    text = f.readlines()
                    self.num_particles = len(text)
                    for i in range(len(text)):
                        self.mass.append((float(text[i].strip('\n')
                                                .split(',')[10]))/1.1)
            else:
                pass


def plotting(num_folders, folders, time_data, mass_data):
    f = plt.figure(figsize=(20, 10))
    ax = f.add_subplot(111)
    for i in range(num_folders):
        if 'an' in folders[i][0:2]:
            ax.plot(time_data[i], mass_data[i], 'k-',
                    label=r'$D^2/{D_0}^2$ for analytic solution')
        else:
            ax.plot(time_data[i], mass_data[i], '--',
                    label=r'$D^2/{D_0}^2$ for $\Delta t =$' +
                    str(timestep_sizes[i]))
        ax.set_xlim(0)
        ax.set_ylim(0, 1)
        plt.xlabel(r't/$\tau_{d0}$')
        plt.ylabel(r'$D^2/{D_0}^2$')
        plt.title('Non-dimensionalised Droplet Diameter')
        plt.legend(loc='upper right')


def output_data(folders, time_data, mass_data,
                processed_data_loc='processed_data//'):
    for i in range(num_folders):
        with open(processed_data_loc + folders[i] +
                  '_verification.txt', 'w') as f:
            time_data[i][::-1]
            f.write('time' + ' ' + 'd2' + ' ' + '\n')
            for j in range(len(time_data[i])):
                f.write(str(time_data[i][j]) + ' ' + str(mass_data[i][j]) +
                        ' ' + '\n')
    time_data[::-1]


folders = [name for name in os.listdir('data')]
num_folders = len(folders)

timestep_sizes = [float(folders[i].split('_')[-3] + '.' +
                        folders[i].split('_')[-2]) for i in range(num_folders)]

time_data = []
mass_data = []

for i in range(num_folders):
    root = 'data//' + folders[i]
    tau = 3.251928
    f = data_processing(root, tau)
    f.sort_files()
    f.read_files()
    time_data.append(f.time_nd)
    mass_data.append(f.mass)

plotting(num_folders, folders, time_data, mass_data)
output_data(folders, time_data, mass_data)
