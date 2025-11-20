#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 15:25:32 2024

@author: philippe.gris@clermont.in2p3.fr
"""
from sn_telmodel.sn_throughputs import Throughputs
from sn_tools.sn_utils import multiproc
import pandas as pd
import numpy as np


class Sigma_zp_meanwave:
    def __init__(self, through_dir, site_name='LSST', pressure=743.,
                 par_names=['airmass', 'pwv', 'ozone', 'beta', 'aerosol'],
                 par_means=[1.2, 4.0, 300., 0.05, 0.05],
                 par_sigmas=[0.01, 0.2, 10., 0.0, 0.001],
                 save_throughputs_dir=''):
        """
        class to estimate sigma_zp and sigma_lambdabar 
        according to atmospheric parameters variation

        Parameters
        ----------
        through_dir : str
            Throughput dir.
        site_name : str, optional
            Site name for observations. The default is 'LSST'.
        pressure : float, optional
            Site pressure. The default is 743..
        par_names : list(str), optional
            parameter list. 
            The default is ['airmass', 'pwv', 'ozone', 'beta', 'aerosol'].
        par_means : list(float), optional
            parameter means. The default is [1.2, 4.0, 300., 0.05, 0.05].
        par_sigmas : list(float), optional
            parameters sigmas. The default is [0.01, 0.2, 10., 0.0, 0.001].

        Returns
        -------
        None.

        """

        self.throughput = Throughputs(
            tel_dir=through_dir,
            site_name=site_name, pressure=pressure)

        self.mean_values = self.get_values(par_names, par_means)
        self.sigma_values = self.get_values(par_names, par_sigmas)

        self.save_throughputs_dir = save_throughputs_dir

        if self.save_throughputs_dir != '':
            from sn_tools.sn_io import checkDir
            checkDir(self.save_throughputs_dir)

    def get_values(self, names, values):
        """
        Method to transform two lists to a dict

        Parameters
        ----------
        names : list(str)
            List of names.
        values : list(float)
            List of values.

        Returns
        -------
        dict
            Resulting dict.

        """

        return dict(zip(names, values))

    def __call__(self, ntrials=1000, nproc=8):
        """
        Main method for data processing

        Parameters
        ----------
        ntrials : int, optional
            number of random parameter choices. The default is 1000.
        nproc : int, optional
            number of procs to use for processing. The default is 8.

        Returns
        -------
        df : pandas df
            Result.

        """

        # get random values
        param_values = self.get_random_values(ntrials)

        idx = param_values['airmass'] >= 1
        idx &= param_values['pwv'] >= 0.
        idx &= param_values['ozone'] >= 0.
        idx &= param_values['aerosol'] >= 0.
        param_values = param_values[idx]
        params = {}

        params['throughput'] = self.throughput
        zp_meanwave = multiproc(param_values, params, self.zp_meanwave, nproc)

        # print('jjj', zp_meanwave.columns)

        cols = zp_meanwave.columns
        fi_vals = {}

        for col in cols:
            means = zp_meanwave[col].mean()
            stds = zp_meanwave[col].std()
            fi_vals['mean_{}'.format(col)] = [means]
            fi_vals['std_{}'.format(col)] = [stds]

        """
        vv = zp_values.mean().to_list()
        cols = zp_values.columns.to_list()
        colsb = list(map(lambda x: 'mean_zp_' + x, cols))
        df = pd.DataFrame([vv], columns=colsb)
        colsc = list(map(lambda x: 'std_zp_' + x, cols))
        df[colsc] = zp_values.std().to_list()
        """
        df = pd.DataFrame.from_dict(fi_vals)

        # add atmospheric parameters
        df = self.concat(df, self.mean_values, 'mean')
        df = self.concat(df, self.sigma_values, 'sigma')

        return df

    def concat(self, dfa, thedict, prefix):
        """
        Method to concat df with dict transformed as df

        Parameters
        ----------
        dfa : pandas df
            original df.
        thedict : dict
            Data to merge.
        prefix : str
            prefix for column names.

        Returns
        -------
        dfa : pandas df
            Resulting merged df.

        """

        dfb = self.make_df(thedict, prefix)
        dfa = pd.concat((dfa, dfb), axis=1)

        return dfa

    def make_df(self, thedict, prefix='mean'):
        """
        Method to create a df from dict with column name change

        Parameters
        ----------
        thedict : dict
            Data to process.
        prefix : str, optional
            prefix to add to col names. The default is 'mean'.

        Returns
        -------
        dfa : pandas df
            Result.

        """

        cols = thedict.keys()
        colsb = list(map(lambda x: '{}_'.format(prefix) + x, cols))
        dfa = pd.DataFrame([thedict.values()], columns=colsb)

        return dfa

    def get_random_values(self, ntrials):
        """
        Method to estimate random values for atmospheric parameters

        Parameters
        ----------
        ntrials : int
            number of random parameter choices.

        Returns
        -------
        res : pandas df
            random values for atmospheric parameters.

        """

        rnd = {}
        for key in self.mean_values.keys():
            rnd[key] = []
        for i in range(ntrials):
            for key, vals in self.mean_values.items():
                sigma = self.sigma_values[key]
                vv = vals+np.random.normal(0., sigma)
                rnd[key].append(vv)

        res = pd.DataFrame.from_dict(rnd)

        return res

    def zp_meanwave(self, data, params, j=0, output_q=None):
        """
        Method to estimate zero points for atmospheric parameters (data)

        Parameters
        ----------
        data : pandas df
            Atmospheric parameters.
        params : dict
            Method parameters.
        j : int, optional
            int for multiprocessing. The default is 0.
        output_q : multiprocessing queue, optional
            Where to store the results. The default is None.

        Returns
        -------
        pandas df
            zp results for each band.

        """

        throughput = params['throughput']
        zp_dict = {}
        bands = list('ugrizy')
        zp_dict = dict(zip(bands, [[], [], [], [], [], []]))
        mean_wave_dict = dict(zip(bands, [[], [], [], [], [], []]))

        df_throughputs = pd.DataFrame()
        for i, row in data.iterrows():
            throughput.reset_data()
            throughput.new_atmosphere(airmass=row['airmass'],
                                      aerosol=row['aerosol'],
                                      pwv=row['pwv'],
                                      ozone=row['ozone'],
                                      beta=row['beta'])
            if self.save_throughputs_dir != '':
                df_throughputs = pd.concat((df_throughputs,
                                            self.get_throughputs(throughput, i, j)))
            throughput.mean_wave()
            for b in 'ugrizy':
                # mean_wave = tel.mean_wavelength[b]
                zpb = throughput.zp(b, exptime=30, nexp=1)
                zp_dict[b].append(zpb)
                mean_wave_dict[b].append(throughput.mean_wavelength[b])

        res_zp = pd.DataFrame.from_dict(zp_dict)
        res_zp.columns = 'zp_' + res_zp.columns
        res_meanwave = pd.DataFrame.from_dict(mean_wave_dict)
        res_meanwave.columns = 'mean_wave_' + res_meanwave.columns
        res = pd.concat((res_zp, res_meanwave), axis=1)

        if self.save_throughputs_dir != '':
            outName = '{}/throughputs_{}.hdf5'.format(
                self.save_throughputs_dir, j)
            df_throughputs.to_hdf(outName, key='throughputs')

        if output_q is not None:
            return output_q.put({j: res})
        else:
            return res

    def get_throughputs(self, throughput, ia, ib):

        df_combi = pd.DataFrame()
        for b in 'ugrizy':
            tt = throughput.throughputs[b]
            dfa = pd.DataFrame(tt.sb, columns=['sb'])
            dfa['wavelen'] = tt.wavelen
            dfa['band'] = b
            df_combi = pd.concat((df_combi, dfa))
            df_combi['combi'] = 'combi_{}_{}'.format(ia, ib)

        return df_combi
