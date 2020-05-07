# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib as mpl
from matplotlib import cm
from scipy import stats
from matplotlib.backends.backend_pgf import FigureCanvasPgf
mpl.backend_bases.register_backend('pdf', FigureCanvasPgf)
import matplotlib.pyplot as plt
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['text.usetex'] = True
mpl.rcParams['pgf.texsystem'] = 'lualatex'
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams.update({'figure.autolayout': True})
from mpl_toolkits.mplot3d import Axes3D


class data_processing():
    def __init__(self, root, tau):
        self.root = root
        self.tau = tau
        self.g = 9.81
        self.time = []
        self.time_nd = []
        self.x_pos = []
        self.y_pos = []
        self.z_pos = []
        self.x_vel = []
        self.y_vel = []
        self.z_vel = []
        self.x_fluid_vel = []
        self.y_fluid_vel = []
        self.z_fluid_vel = []
        self.ang_bar = []
        self.x_vel_bar = []
        self.y_vel_bar = []
        self.z_vel_bar = []
        self.T_d_prime = []
        self.temp_mean_data = []
        self.temp_std_data = []
        self.mass_mean_data = []
        self.mass_std_data = []
        self.T_d = []
        self.m_d = []
        self.kd = []
        self.kd_vals = []
        self.T_d_mean = 0
        self.num_particles = 0

        self.temp_dx = []
        self.temp_dy = []
        self.temp_dz = []
        self.temp_x_vals = []
        self.temp_y_vals = []
        self.temp_z_vals = []
        self.temp_rgba = []

        self.mass_dx = []
        self.mass_dy = []
        self.mass_dz = []
        self.mass_x_vals = []
        self.mass_y_vals = []
        self.mass_z_vals = []
        self.mass_rgba = []

        self.temp_data = []
        self.mass_data = []

        self.temp_pdf_data = []
        self.mass_pdf_data = []
        self.temp_y_vals = []
        self.mass_y_vals = []

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
        #self.files = self.files[::10]
        self.files = self.files[0:801] # 100

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
                        self.T_d.append(float(text[i].strip('\n')
                                              .split(',')[9]))
                        self.m_d.append(float(text[i].strip('\n')
                                              .split(',')[10]))
            else:
                pass

    def temp_data_gen(self):
        i = 0
        while i < len(self.T_d):
            self.T_d_mean = np.mean(self.T_d[i])
            self.T_d_prime.append([self.T_d[i] - self.T_d_mean
                                   for i in range(i, i+self.num_particles)])
            i += self.num_particles

        for i in range(len(self.T_d_prime)):
            self.temp_mean_data.append(np.mean(self.T_d_prime[i]))
            self.temp_std_data.append(np.std(self.T_d_prime[i]))

    def mass_data_gen(self):
        self.m_d_0 = self.m_d[0]
        print(self.m_d[0])

        i = 0
        while i < len(self.m_d):
            self.mass_data.append([self.m_d[i]/self.m_d_0 for i in
                                   range(i, i+self.num_particles)])
            i += self.num_particles

        for i in range(len(self.T_d_prime)):
            self.mass_mean_data.append(np.mean(self.mass_data[i]))
            self.mass_std_data.append(np.std(self.mass_data[i]))

    def temp_hist_generator(self):
        for i in range(len(self.T_d_prime)):
            self.hist, self.xedges = np.histogram(self.T_d_prime[i],
                                                  bins='auto', density=True)

            self.temp_dx.append(self.xedges[1] - self.xedges[0])
            self.temp_dy.append(self.time[1] - self.time[0])
            self.temp_dz.append(self.hist.flatten())
            self.temp_x_vals.append(np.meshgrid(self.xedges[:-1] +
                                                self.xedges[1:])[0]/2)
            self.temp_y_vals.append([self.time_nd[i]]*len(self.temp_x_vals[i]))
            self.temp_z_vals.append(np.zeros_like(self.temp_x_vals[i]))

            self.cmap = cm.get_cmap('jet')
            self.max_height = np.max(self.temp_dz[i])
            self.min_height = np.min(self.temp_dz[i])
            # scale each z to [0,1], and get their rgb values
            self.temp_rgba.append([self.cmap((k-self.min_height) /
                                   self.max_height) for k in self.temp_dz[i]])

    def mass_hist_generator(self):
        for i in range(len(self.mass_data)):
            self.hist, self.xedges = np.histogram(self.mass_data[i],
                                                  bins='auto', density=True)

            self.mass_dx.append(self.xedges[1] - self.xedges[0])
            self.mass_dy.append(self.time[1] - self.time[0])
            self.mass_dz.append(self.hist.flatten())
            self.mass_x_vals.append(np.meshgrid(self.xedges[:-1] +
                                                self.xedges[1:])[0]/2)
            self.mass_y_vals.append([self.time_nd[i]] * len(self.x_vals[i]))
            self.mass_z_vals.append(np.zeros_like(self.x_vals[i]))

            self.cmap = cm.get_cmap('jet')
            self.max_height = np.max(self.dz[i])
            self.min_height = np.min(self.dz[i])
            # scale each z to [0,1], and get their rgb values
            self.mass_rgba.append([self.cmap((k-self.min_height) /
                                   self.max_height) for k in self.dz[i]])

    def temp_pdf(self):
        for i in range(len(self.T_d_prime)):
            self.temp_pdf_data.append(stats.norm.pdf(self.T_d_prime[i],
                                                     self.temp_mean_data[i],
                                                     self.temp_std_data[i]))

            self.temp_y_vals.append([self.time_nd[i]] * len(self.T_d_prime[i]))

    def mass_pdf(self):
        for i in range(len(self.mass_data)):
            self.mass_pdf_data.append(stats.norm.pdf(self.mass_data[i],
                                                     self.mass_mean_data[i],
                                                     self.mass_std_data[i]))

            self.mass_y_vals.append([self.time_nd[i]] * len(self.mass_data[i]))


def plot_temp_pdfs(T_d_prime, temp_y_vals, temp_pdf_data, figsave=False):
    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(111, projection='3d')
    for i in range(0, 50):
        ax.scatter(T_d_prime[i][::10], temp_y_vals[i][::10],
                   temp_pdf_data[i][::10], '.', s=1,
                   c=temp_pdf_data[i][::10], cmap='jet')
    ax.set_xlabel(r'$T^\prime$', fontsize=35, labelpad=40)
    ax.set_ylabel(r'$t/\tau_{d0}$', fontsize=35, labelpad=40, ha='right')
    ax.set_zlabel(r'Density', fontsize=35, labelpad=40)
    ax.tick_params(labelsize=25)
    ax.xaxis.set_rotate_label(False)
    ax.yaxis.set_rotate_label(False)
    ax.zaxis.set_rotate_label(False)
    plt.xlim(-0.03, 0.03)
    plt.ylim(0)
    ax.set_zlim(0)
    ax.view_init(elev=10, azim=-45)
    if figsave is True:
        plt.savefig('Temperature_PDFs.png', dpi=300, bbox_inches='tight',
                    pad_inches=1)


def plot_mass_pdfs(mass_data, mass_y_vals, mass_pdf_data, figsave=False):
    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(111, projection='3d')
    for i in range(0, 50):
        ax.scatter(mass_data[i][::10], mass_y_vals[i][::10],
                   mass_pdf_data[i][::10], '.', s=1,
                   c=mass_pdf_data[i][::10], cmap='jet')
    ax.set_xlabel(r'$D^2/{D_0}^2$', fontsize=35, labelpad=40)
    ax.set_ylabel(r'$t/\tau_{d0}$', fontsize=35, labelpad=40)
    ax.set_zlabel(r'Density', fontsize=35, labelpad=50)
    ax.tick_params(labelsize=25)
    ax.xaxis.set_rotate_label(False)
    ax.yaxis.set_rotate_label(False)
    ax.zaxis.set_rotate_label(False)
    plt.xlim(0.95, 1)
    plt.ylim(0)
    ax.set_zlim(0)
    ax.view_init(elev=15, azim=-100)
    if figsave is True:
        plt.savefig('Mass_PDFs.png', dpi=300, bbox_inches='tight',
                    pad_inches=1)


def plot_temp_mass_pdf(mass_data, T_d_prime, timestep_no=800, figsave=False):
    f1 = plt.figure(figsize=(12, 12))
    ax = f1.add_subplot(111)
    ax.plot(mass_data[timestep_no], T_d_prime[timestep_no], '.')
    ax.tick_params(labelsize=25)
    ax.set_xlabel(r'$D^2/{D_0}^2$', fontsize=35, labelpad=40)
    ax.set_ylabel(r'$T^\prime$', fontsize=35, labelpad=40)
    if figsave is True:
        plt.savefig('Temp_Mass_PDFs_timestep_' + str(timestep_no) +
                    '.png', dpi=300, bbox_inches='tight',
                    pad_inches=1)


def plot_deviation_range(time_data, var_data, T_d_prime, figsave=True):
    stks_nums = ['0_1', '1_0', '10']
    for i in range(len(time_data)):
        f1 = plt.figure(figsize=(10, 10))
        ax = f1.add_subplot(111)
        ax.plot(time_data[i][0:50], var_data[i][0:50], '--')
        ax.tick_params(labelsize=25)
        ax.set_xlabel(r'$t/\tau_{d0}$', fontsize=35, labelpad=40)
        ax.set_ylabel(r'$T^\prime~Range$', fontsize=35, labelpad=40)
        if figsave is True:
            plt.savefig('Temp_Deviation_' + stks_nums[i] + '.png', dpi=300,
                        bbox_inches='tight', pad_inches=1)


def plot_histograms():
    root = 'D:\Andrew\Documents\Andrew University\Part III\Individual Project\HMOranges\cmake-build-debug\TGV_PERIODIC_1'#'data//' + folders[i]
    tau = (2000*0.001048**2)/(18*0.000018735)
    f = data_processing(root, tau)
    f.sort_files()
    f.read_files()
    f.temp_data_gen()
    f.temp_pdf()
    f.mass_data_gen()
    f.mass_pdf()
    f.temp_hist_generator()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(1, len(f.mass_data)-1):
        ax.bar3d(f.x_vals[i], f.y_vals[i], f.z_vals[i], f.dx[i], f.dy[i],
                 f.dz[i], color=f.rgba[i], zsort='average', shade=True)
    plt.title('Histogram of Temporal Temperature Data')
    plt.xlabel(r'$T\prime$')
    plt.ylabel(r'$t/\tau_d$')


def plot_PDFs():
    time_data = []
    #for i in range(num_folders):
    root = 'D:\Andrew\Documents\Andrew University\Part III\Individual Project\HMOranges\cmake-build-debug\TGV_PERIODIC_1'#'data//' + folders[i]
    tau = (2000*0.001048**2)/(18*0.000018735)
    f = data_processing(root, tau)
    f.sort_files()
    f.read_files()
    f.temp_data_gen()
    f.temp_pdf()
    f.mass_data_gen()
    f.mass_pdf()
    time_data.append(f.time_nd)

    plot_mass_pdfs(f.mass_data, f.mass_y_vals, f.mass_pdf_data, False)
    plot_temp_pdfs(f.T_d_prime, f.temp_y_vals, f.temp_pdf_data, False)
    plot_temp_mass_pdf(f.mass_data, f. T_d_prime)


def plot_deviation():
    folders = [name for name in os.listdir('stokes_numbers\data')]
    mu_G = [2.2804*10**(-5), 2.2804*10**(-6), 2.2804*10**(-7)]
    time_data = []
    var_data = []
    for i in range(len(folders)):
        root = 'stokes_numbers\data//' + folders[i]
        tau = (997*0.001048**2)/(18*mu_G[i])
        f = data_processing(root, tau)
        f.sort_files()
        f.read_files()
        f.temp_data_gen()
        time_data.append(f.time_nd)
        var_data.append([np.max(f.T_d_prime[i]) - np.min(f.T_d_prime[i])
                         for i in range(len(f.T_d_prime))])
    plot_deviation_range(time_data, var_data, f.T_d_prime)


#plot_PDFs()
plot_deviation()

#x = [item for sublist in f.T_d_prime[0:10] for item in sublist]
#y = [item for sublist in f.y_vals[0:10] for item in sublist]
#z = [item for sublist in f.temp_pdf_data[0:10] for item in sublist]
#
#x = x[::100]
#y = y[::100]
#z = z[::100]

#x = [f.temp_pdf_data[i][j] for i in range(len(f.time)) for j in range(len(f.temp_pdf_data[i]))]
#x = []
#z = []
#j = 0
#while j < (len(f.T_d_prime[0])):
#    x.append([f.T_d_prime[i][j] for i in range(len(f.T_d_prime))])
#    z.append([f.temp_pdf_data[i][j] for i in range(len(f.time))])
#    j += 10

#for i in range(0,1000):
#    ax.plot(x[i], f.time, z[i], '-', color='white')#cmap[i])



#plt.savefig('test.pdf')
#for ii in range(0,360,1):
#    ax.view_init(elev=10., azim=ii)
#    plt.draw()
#    plt.pause(.001)





#ax.bar(xpos, dz, width=(xpos[1]-xpos[0]), align='center')

#f1 = plt.figure()
#ax1 = f1.add_subplot(111)
#ax1.hist(f.T_d_prime[190], bins='auto', density=0)
#plt.show
