#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 10:19:42 2024

@author: philippe.gris@clermont.in2p3.fr
"""
from rubin_sim.phot_utils import Bandpass
import os


class Telescope:
    def __init__(self, tel_dir='throughputs_1.9/baseline',
                 tel_optical_files=['detector.dat', 'lens1.dat',
                                    'lens2.dat', 'lens3.dat',
                                    'm1.dat', 'm2.dat', 'm3.dat'],
                 tel_filter_files=['filter_u.dat', 'filter_g.dat',
                                   'filter_r.dat', 'filter_i.dat',
                                   'filter_z.dat', 'filter_y.dat'],
                 tel_wave_min=300.,
                 tel_wave_max=1150.,
                 filter_colors=dict(zip(['u', 'g', 'r', 'i', 'z', 'y'],
                                        ['b', 'c', 'g', 'y', 'r', 'm']))):
        """
        Telescope class

        Parameters
        ----------
        tel_dir : str, optional
            input files directory. The default is 'throughputs_1.9/baseline'.
        tel_optical_files : list(str), optional
            List of optical components. The default is
                                        ['detector.dat', 'lens1.dat',
                                         'lens2.dat', 'lens3.dat',
                                         'm1.dat', 'm2.dat', 'm3.dat'].
        tel_filter_files : list(str), optional
            List of filter files. The default is
                                        ['filter_u.dat', 'filter_g.dat',
                                         'filter_r.dat', 'filter_i.dat',
                                        'filter_z.dat', 'filter_y.dat'].
        tel_wave_min : float, optional
            Min wave length. The default is 300..
        tel_wave_max : float, optional
            Max wavelength. The default is 1150..
       filter_colors : list(str), optional
            List of filter colors (for display). The default is 
                dict(zip(['u', 'g', 'r', 'i', 'z', 'y'],                                        
                         ['b', 'c', 'g', 'y', 'r', 'm'])).

        Returns
        -------
        None.

        """

        self.tel_dir = tel_dir
        self.tel_optical_files = tel_optical_files
        self.tel_filter_files = tel_filter_files
        self.tel_wave_min = tel_wave_min
        self.tel_wave_max = tel_wave_max
        self.filter_colors = filter_colors

        self.filter_list = self.get_filter_list()

        self.load_components()
        self.tel_trans = {}

        self.load_system()

    def get_filter_list(self):
        """
        Method to estimate the filter list from the tel_filter_files

        Returns
        -------
        flist : str
            Filter list.

        """

        flist = ''
        for fi in self.tel_filter_files:
            f = fi.split('.dat')[0].split('_')[-1]
            flist += f

        return flist

    def load_components(self):
        """
        Method to load all components (detector, lenses, filters, mirrors)

        Returns
        -------
        None.

        """

        self.optical = {}

        for fi in self.tel_optical_files:
            vv = Bandpass()
            path_vv = os.path.join(self.tel_dir, fi)
            print('reading', path_vv)
            vv.read_throughput(path_vv)
            nname = fi.split('.dat')[0]
            self.optical[nname] = vv

        self.filter = {}
        for fi in self.tel_filter_files:
            vv = Bandpass()
            path_vv = os.path.join(self.tel_dir, fi)
            vv.read_throughput(path_vv)
            f = fi.split('.dat')[0].split('_')[-1]
            self.filter[f] = vv

    def load_system(self):
        """
        Method to estimate the total telescope transmission (optics+filters)
        The resul is in the dict self.tel_trans[b] where b is the band

        Returns
        -------
        None.

        """

        for fi in self.filter_list:
            f = fi.split('.dat')[0].split('_')[-1]
            self.tel_trans[f] = Bandpass()

            index = [i for i, x in enumerate(
                self.tel_filter_files) if f+'.dat' in x]
            telfiles = self.tel_optical_files + \
                [self.tel_filter_files[index[0]]]

            print(f, telfiles)
            self.tel_trans[f].read_throughput_list(telfiles,
                                                   root_dir=self.tel_dir,
                                                   wavelen_min=self.tel_wave_min,
                                                   wavelen_max=self.tel_wave_max)

    def plot_components(self):
        """
        Methos to plot telescope components

        Returns
        -------
        None.

        """

        # optical components
        for key, vals in self.optical.items():
            self.plot_component(key, vals, figtit='Telescope {}'.format(key))

        import matplotlib.pyplot as plt
        # filters
        fig, ax = plt.subplots(figsize=(12, 8))
        # without optics
        for key, vals in self.filter.items():
            self.plot_component(key, vals, fig=fig, ax=ax,
                                color=self.filter_colors[key],
                                label='{} band'.format(key))

        # with optics
        for key, vals in self.tel_trans.items():
            self.plot_component(key, vals, fig=fig, ax=ax,
                                color=self.filter_colors[key],
                                linestyle='dashed',
                                label='{} band+optics'.format(key))

    def plot_component(self, compo, trans, fig=None, ax=None,
                       color='black', linestyle='solid', figtit='',
                       label=''):
        """
        Method to plot a single component

        Parameters
        ----------
        compo : str
            Component to plot.
        trans : transmission (.wl and .sb)
            Transmission.
        fig : matplotlib figure, optional
            Figure for the plot. The default is None.
        ax : matplotlib axis, optional
            Axis for the plot. The default is None.
        color : str, optional
            Color for the plot. The default is 'black'.
        linestyle : str, optional
            Linestyle for the plot. The default is 'solid'.
        figtit: str, optional
            Figure title. The default is ''
        label: str, optional
            label for the plot. The default is ''

        Returns
        -------
        None.

        """

        import matplotlib.pyplot as plt
        if fig is None:
            fig, ax = plt.subplots(figsize=(12, 8))

        ax.plot(trans.wavelen, trans.sb, linestyle=linestyle,
                color=color, label=label)
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Sb (0-1)')
        ax.set_title(figtit)
        ax.set_ylim([0.0, 1.])
        ax.grid(visible=True)
        ax.legend()
