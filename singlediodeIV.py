#!/usr/bin/env/ python3

import scipy.special.lambertw as W
import numpy as np
import matplotlib.pyplot as plt

# note: the docstrings are indented so that they can be easily folded

def generateIV(*user_parameters):
    """
        Generates current-voltage data using the Shockley diode equation.

        The default parameters will return dark current-voltage data for a diode
        with n = 1, T = 300 K, I0 = 1 nA, IL = 0 A, Rs = 10 ohm, and Rsh = 1 Mohm

        The area defaults to unity; if an area is specified, the function returns
        current density-voltage instead.

        Parameters
        ----------
        *user_parameters: list
            must be in format [n, T, I0, Iph, Rs, Rsh]

        n : numeric
            diode ideality factor

        T : numeric
            absolute temperature [K]

        I0 : numeric
            diode reverse saturation current [A]

        Iph : numeric
            photo-generated current [A]

        Rs : numeric
            diode series resitance [ohms]

        Rsh : numeric
            diode shunt resistance [ohms]


        Returns
        -------
        list
            [[voltages [V]], [currents [A]]
    """
    default = [1, 300, 1E-9, 0.0, 10, 1E6]
    # if the user submits a parameter list, parameters will update
    if user_parameters:
        parameters = []
        for i in user_parameters:
            parameters += i
        print('\n'+'Collected user input parameters:'+'\n'+
              '{}'.format(parameters))
    else:
        parameters = default
        print('\n'+'Data generated using default parameters:'+'\n'+
              '{}'.format(parameters))

    # set all variables to their corresponding parameters
    n, T, I0, Iph, Rs, Rsh = parameters

    # generate a default voltage range for now
    step_size = 0.01
    voltages = list(np.arange(-2.0, 1.0, step_size))

    # initialize the current list
    currents = []

    # store the thermal voltage [V]
    Vth = 8.617332E-5*T

    # calculate current for each voltage
    for V in voltages:
        z = (Rs*I0)/(n*Vth)*np.exp((Rs*(Iph+I0)+V)/(n*Vth*(1+Rs/Rsh)))
        I = ((Iph + I0) - V/Rsh)/(1 +Rs/Rsh) - (n*Vth)/(Rs) * W(z, k=0)
        currents.append(-I.real)
        V += step_size

    IV_data = [voltages, currents]
    return IV_data


def make_fig(*data, **kwargs):
    """
        Outputs a matplotlib.figure.Figure object.

        *args
        -----
        data : list
            Must be a list containing two lists (for x and y)

        *kwargs
        -------
        title : string
            Creates the figure with a title

        xlabel : string
            Labels the x-axis

        ylabel : string
            Labels the y-axis
    """
    # initialize the figure with an empty subplot
    fig, ax = plt.subplots(1, 1, figsize=(5,5))

    # handling *args
    if data:
        xdata = data[0][0]
        ydata = data[0][1]
        ax.plot(xdata, ydata)
    else:
        ax.plot()

    # handling *kwargs
    for name, value in kwargs.items():
        if name == 'title':
            ax.set_title(value)
        elif name == 'xlabel':
            ax.set_xlabel(value)
        elif name == 'ylabel':
            ax.set_ylabel(value)

    return fig


def make_semilog_fig(*data, **kwargs):
    """
        Docstring goes here.
    """
    # initialize the figure with an empty subplot
    fig, ax = plt.subplots(1, 1, figsize=(5,5))

    # handling *args
    if data:
        xdata = data[0][0]
        ydata = [abs(val) for val in data[0][1]]
        ax.plot(xdata, ydata)
    else:
        ax.plot()

    # set the y-axis of the subplot to log_10 scale
    ax.set_yscale('log')

    # handling *kwargs
    for name, value in kwargs.items():
        if name == 'title':
            ax.set_title(value)
        elif name == 'xlabel':
            ax.set_xlabel(value)
        elif name == 'ylabel':
            ax.set_ylabel(value)

    return fig


# this one is not working
def setZoom(figure, xmin=-0.1, xmax=1.0, ymin=-0.03, ymax=0.01):
    plt.axis([xmin, xmax, ymin, ymax])
    return figure


# this should be re-written to take two Figure objects as arguments
def plotIV_inline(data, title='default title'):

    abs_current = [abs(val) for val in data[1]]

    figure, (ax1, ax2) = plt.subplots(1,2)
    figure.suptitle(title)
    figure.set_figwidth(15)
    figure.set_figheight(5)


    ax1.plot(data[0], abs_current)
    ax1.set(title = 'abs(I) vs. V [semilog]',
            xlabel = 'voltage / [V]',
            ylabel = 'current / [A]',
            yscale = 'log',
            ylim = ([1.0E-9, 10]),
            )


    ax2.plot(data[0], data[1])
    ax2.set(title = 'linear I vs. V [linear]',
            xlabel='voltage / [V]',
            ylabel='current / [A]',
            )

    ax2.axhline(0, color='black', linewidth=1.0)
    ax2.axvline(0, color='black', linewidth=1.0)

    return figure



