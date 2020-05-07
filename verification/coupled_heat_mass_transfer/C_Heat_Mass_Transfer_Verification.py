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
    def __init__(self, root):
        self.root = root
        self.time = []
        self.temp = []
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
                    text = f.readlines()
                    self.num_particles = len(text)
                    for i in range(len(text)):
                        self.temp.append(float(text[i].strip('\n')
                                               .split(',')[9]))
                        self.mass.append(float(text[i].strip('\n')
                                               .split(',')[10]))
            else:
                pass

        self.time = np.asarray(self.time, dtype=np.float64)
        self.temp = np.asarray(self.temp, dtype=np.float64)
        self.mass = np.asarray(self.mass, dtype=np.float64)

        self.time = self.time[self.mass > 0]
        self.temp = self.temp[self.mass > 0]
        self.mass = self.mass[self.mass > 0]


def plotting(num_folders, folders, time_data, temp_data, mass_data):
    for i in range(num_folders):
        fig, (ax1, ax2) = plt.subplots(2, figsize=(20, 10))
        ax1.plot(time_data[i], mass_data[i], '--')
        ax1.set_xlim(0)
        ax1.set_ylim(0)
        ax1.set_xlabel(r'$t$ ($s$)')
        ax1.set_ylabel(r'$D^2$ ($mm^2$)')
        ax1.set_title('Diameter Evolution of Evaporating ' +
                      folders[i].split('_')[-1].title() + ' Droplet')

        ax2.plot(time_data[i], temp_data[i], '--')
        ax2.set_xlim(0)
        ax2.set_ylim(temp_data[i][0])
        ax2.set_xlabel(r'$t$ ($s$)')
        ax2.set_ylabel(r'$T_d$ ($K$)')
        ax2.set_title('Temperature Evolution of Evaporating ' +
                      folders[i].split('_')[-1].title() + ' Droplet')


def output_data(folders, time_data, temp_data, mass_data,
                processed_data_loc='processed_data//'):
    for i in range(num_folders):
        with open(processed_data_loc + folders[i] +
                  '_verification.txt', 'w') as f:
            time_data[i][::-1]
            f.write('time' + ' ' + 'T_d' + ' ' + 'd2' + ' ' + '\n')
            for j in range(len(time_data[i])):
                f.write(str(time_data[i][j]) + ' ' + str(temp_data[i][j]) +
                        ' ' + str(mass_data[i][j]) + ' ' + '\n')
    time_data[::-1]


folders = [name for name in os.listdir('data')]
num_folders = len(folders)

time_data = []
temp_data = []
mass_data = []

for i in range(num_folders):
    root = 'data//' + folders[i]
    f = data_processing(root)
    f.sort_files()
    f.read_files()
    time_data.append(f.time)
    temp_data.append(f.temp)
    mass_data.append(f.mass)

plotting(num_folders, folders, time_data, temp_data, mass_data)
output_data(folders, time_data, temp_data, mass_data)
