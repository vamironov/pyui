#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 08:38:06 2019

@author: Vladimir Mironov
"""

import argparse
import numpy as np
import pandas as pd
import scipy.integrate
import scipy.interpolate

KB = 0.0019872041  # kcal/mol


def us_umbrella_integration(data, win_centers, win_k,
                            temperature, xmin, xmax, nbins, verbose=False):
    """
    Run umbrella integration analysis

    Parameters
    ----------
    data : pandas.DataFrame
        data frame containing statistical data of the reaction coordinate value
        for every US window in separate columns

    win_centers : numpy.array(float64)
        centers of umbrella potential in every window

    win_k : list of float
        force constants of umbrella potential in every window

    temperature : float
        simulation temperature in K

    xmin : float
        minimum value of reaction coordinate in histogram

    xmax : float
        maximum value of reaction coordinate in histogram

    nbins : int
        number of histogram bins

    verbose : bool, optional
        whether to print windows statistics, default False
    """

    beta_inv = (KB*temperature)
    mean = np.array(data.mean(numeric_only=True))
    var = np.array(data.var(ddof=0, numeric_only=True))
    nval = np.array(data.count(numeric_only=True))

    xstep = (xmax-xmin)/nbins
    xbins = np.linspace(xmin, xmax, nbins+1)[:-1] + xstep/2  # centers of bins

    if verbose:
        print('# Windows stats:')
        print(('# {:>4}   {:>12}{:>16}{:>12}{:>10}').format('ID', 'mean', 'var.', 'f.c.', 'n.pts.'))
        for i, (m, v, n, fc) in enumerate(zip(mean, var, nval, win_k)):
            print(f'# [{i:3}]: {m:12.5f}{v:16.5e}{fc:12.5f}{n:10}')

    da_bins = list()

    for x in xbins:
        p_bin = nval/np.sqrt(2*np.pi*var) * np.exp(-0.5*(x-mean)**2/var)
        pnorm_bin = sum(p_bin)
        da_bin = beta_inv * (x-mean) / var - win_k * (x-win_centers)
        da_bins.append(sum(p_bin * da_bin)/pnorm_bin)

    da = scipy.interpolate.CubicSpline(xbins, da_bins)
    a = da.antiderivative()

    return xbins, a, da


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Umbrella integration ' +
                                     'analysis of umbrella ' +
                                     'sampling simulation')

    parser.add_argument('-c', '--config', required=True,
                        type=argparse.FileType('r'),
                        help='name of configuration file')

    parser.add_argument('-m', '--xmin', type=str,
                        help='minimum value of interest, float', default='auto')

    parser.add_argument('-M', '--xmax', type=str,
                        help='maximum value of interest, float', default='auto')

    parser.add_argument('-n', '--nbins', required=True, type=int,
                        help='number of bins')

    parser.add_argument('-t', '--temp', required=True, type=float,
                        help='simulation temperature')

    parser.add_argument('-u', '--units', choices=['kj', 'kcal'], type=str,
                        help='force constant units of energy', default='kcal')

    parser.add_argument('-v', '--verbose', action='store_true',
                                help='verbose output')

    parser.add_argument('-d', '--data_hist', action='store_true',
                        help='print histogram for each window')

    args = parser.parse_args()

    if args.units == 'kj':
        KB = KB*4.184

    opt_pack = np.genfromtxt(args.config.name,
                             dtype=['U255', np.float64, np.float64])

    files, wlist, flist = map(list, zip(*opt_pack))

    dlist = [pd.read_table(f, names=['t', 'v'], sep=r'\s+',
                           usecols=['v']) for f in files]
    df = pd.concat(dlist, ignore_index=True, axis=1)

    data_min = df.min().min()
    data_max = df.max().max()

    xmin = data_min if args.xmin == 'auto' else float(args.xmin)
    xmax = data_max if args.xmax == 'auto' else float(args.xmax)

    if args.data_hist:
        bins = np.linspace(xmin, xmax, args.nbins+1)
        hists = pd.DataFrame()
        for s in df:
            h, be = np.histogram(df[s], bins=bins)
            hists[s] = pd.Series(1.0*h/h.sum(), index=be[:-1])

        print('# Data distribution for umbrella sampling windows')
        print(hists.to_string())
        print()


    if args.verbose:
        print('# US parameters:')
        print(f'# Xmin  = {xmin:+.5f}')
        print(f'# Xmax  = {xmax:+.5f}')
        print(f'# Nbins = {args.nbins}')
        print(f'# Temp  = {args.temp}')
        print(f'# Units = {args.units}')


    xb, eb, _ = us_umbrella_integration(df, wlist, flist, args.temp,
                                         xmin, xmax, args.nbins, args.verbose)

    E_STR = 'E, '+args.units+'/mol'
    print('X'.rjust(16) + E_STR.rjust(20))

    E_MIN = min(eb(xb))

    for x in xb:
        print("%16.4f%20.8e" % (x, eb(x)-E_MIN))
