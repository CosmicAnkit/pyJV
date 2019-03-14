#!/usr/bin/env python

import scipy.special.lambertw as W
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# example format for parameter input for SingleDiode class
testing_params = {'n': 1.66,
                  'T': 200,
                  'I_0': 3.32E-9,
                  'I_ph': 0.0,
                  'R_s': 9.23,
                  'R_sh': 6.36E3}


class SingleDiode(object):

    def __init__(self, parameters):
        self.ideality_factor = parameters.get('n', 1.0)
        self.temperature = parameters.get('T', 300.0)
        self.saturation_current = parameters.get('I_0', 1.0E-9)
        self.photocurrent = parameters.get('I_ph', 0.0)
        self.resistance_series = parameters.get('R_s', 10.0)
        self.resistance_shunt = parameters.get('R_sh', 1.0E6)
        self.area = parameters.get('area', 1.0)

    def __repr__(self):
        return 'Diode with parameters:\n {}'.format(self.parameters())

    # class methods

    def parameters(self):
        return('   n = {}\n'.format(self.ideality_factor) +
               '   T = {} K\n'.format(self.temperature) +
               ' I_0 = {} amp\n'.format(self.saturation_current) +
               'I_ph = {} amp\n'.format(self.photocurrent) +
               ' R_s = {} ohm\n'.format(self.resistance_series) +
               'R_sh = {} ohm\n'.format(self.resistance_shunt)+
               'area = {} cm^2\n'.format(self.area))

    def generateIV(self, Vmin=-2.0, Vmax=1.0, step=0.01):
        '''
            Generates voltage and current data using the Shockley diode
             equation.
            Diode parameters are obtained from the SingleDiode object.
            Voltage minimum, maximum, and step size are optional inputs.
        '''

        # initialize an empty list for current and generate the voltages
        self.currents = []
        self.voltages = list(np.arange(Vmin, Vmax, step))

        # set variables from object parameters
        n   = self.ideality_factor
        T   = self.temperature
        I0  = self.saturation_current
        Iph = self.photocurrent
        Rs  = self.resistance_series
        Rsh = self.resistance_shunt
        A   = self.area

        # thermal voltage is k_B*T
        Vth = 8.617332E-5 * T

        for V in self.voltages:
            z = (Rs*I0)/(n*Vth)*np.exp((Rs*(Iph+I0)+V)/(n*Vth*(1+Rs/Rsh)))
            I = ((Iph + I0) - V/Rsh)/(1 +Rs/Rsh) - (n*Vth)/(Rs) * W(z, k=0)
            self.currents.append(-I.real)
            V += step

        self.IV_data = [self.voltages, self.currents]
        self.JV_data = [self.voltages, [I/A for I in self.currents]]
        return self.IV_data

    # Need to make this function extendable (not class method)
    def plotIV(self):
        lin_fig = plt.figure(num='Linear Current vs. Voltage')
        plt.plot(self.voltages, self.currents)
        plt.text(-1.5, 0.01, self.parameters(), fontname='monospace')
        plt.ylabel('current / [A]')
        plt.xlabel('voltage / [V]')
        plt.axhline(color='black', linewidth=0.5)
        plt.axvline(color='black', linewidth=0.5)

        log_fig = plt.figure(num='Semi-log Current vs. Voltage')
        plt.plot(self.voltages, [abs(I) for I in self.currents])
        plt.yscale('log')
        plt.ylabel('|current| / [A]')
        plt.xlabel('voltage / [V]')
        # ax = plt.axes(x)
        # ax.set_major_formatter(ScalarFormatter())
        # ax.grid()
        plt.axvline(color='black', linewidth=0.5)

        return lin_fig, log_fig
