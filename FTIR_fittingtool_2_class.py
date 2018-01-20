import os
# import time
import math
from tkinter import *
from tkinter import filedialog, messagebox
from tkinter.ttk import Progressbar
import csv
import matplotlib
matplotlib.use("TkAgg")
# import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import numpy as np
from scipy.optimize import curve_fit
import threading
import queue
import cross_platform_config
from sys import platform as _platform

__version__ = '2.40'


class FIT_FTIR:
    def __init__(self, wavenumbers, transmissions, subd, layertype_list, entry_x_list, entry_d_list,
                 checklayer_list, scalefactor, angle, CdTe_offset, HgTe_offset, subtype, fittype, listbox, progress_var, wn_beingcalculated):
        self.wns = wavenumbers
        self.trans = transmissions
        self.subd = subd
        self.layertype_list, self.entry_x_list, self.entry_d_list, self.checklayer_list \
            = layertype_list, entry_x_list, entry_d_list, checklayer_list
        self.scalefactor = scalefactor
        self.angle = angle
        self.CdTe_offset = CdTe_offset
        self.HgTe_offset = HgTe_offset
        self.subtype = subtype
        self.fittingtype = fittype
        self.listbox = listbox
        self.progress_var = progress_var
        self.wn_beingcalculated = wn_beingcalculated

        self.n_list = []
        self.k_list = []
        self.etas_list = []
        self.etap_list = []
        self.deltas_list = []
        self.deltap_list = []
        self.matrixs_list = []
        self.matrixp_list = []

        self.transs = []
        self.peakvalues = []

        self.alphaMCT = 0
        self.alphaCT = 0

        self.basek = 0
        self.ab = 0
        self.reflections = []
        self.absorptions = []

        self.load_n_file()

        if self.fittingtype == 1:
            self.fit_curve()
        elif self.fittingtype == 2:
            self.show_fringes()
        elif self.fittingtype == 0:
            pass

    def fit_curve(self):
        for wn in self.wns:
            self.lamda = 10000 / float(wn)
            self.C1 = 0.5694 * math.pow(self.x, -1.5355)
            if self.lamda < self.C1 * 1.17:
                self.lamda = self.C1 * 1.17
            self.E = 4.13566743 * 3 / 10 / self.lamda

            self.R1 = self.cal_R(-1, self.x, self.lamda)
            self.R2 = self.cal_R(self.x, 1, self.lamda)
            self.R3 = self.cal_R(1, -1, self.lamda)

            self.cal_alpha_all()

            self.transs.append(float(self.cal_trans() * 100))

        self.addlog('Fitting curve complete.')

    def show_fringes(self):
        self.peakvalues = []
        self.reflections = []
        self.absorptions = []
        for wn in self.wns:
            self.lamda = 10000 / float(wn)
            self.E = 4.13566743 * 3 / 10 / self.lamda
            self.peakvalues.append(self.cal_fringes_single(self.lamda)[0])
            self.reflections.append(self.cal_fringes_single(self.lamda)[1])
            self.absorptions.append(self.cal_fringes_single(self.lamda)[2])

    def cal_fringes_single(self, lamda):
        self.n_list = []
        self.k_list = []
        self.etas_list = []
        self.etap_list = []
        self.deltas_list = []
        self.deltap_list = []
        self.matrixs_list = []
        self.matrixp_list = []

        self.allresult = []

        self.eta0s = np.cos(self.angle)
        self.eta0p = 1 / np.cos(self.angle)

        for i in range(0, len(self.layertype_list)):
            if self.layertype_list[i] == "CdTe":
                self.cal_initialpara(1)
            elif self.layertype_list[i] == "MCT" or self.layertype_list[i] == "SL":
                self.cal_initialpara(self.entry_x_list[i])
            n = self.cal_n(lamda, self.layertype_list[i])
            k = self.cal_k(lamda, self.layertype_list[i])
            self.n_list.append(n)
            self.k_list.append(k)
            etas = np.sqrt((n - 1j * k) * (n - 1j * k) - np.sin(self.angle) * np.sin(self.angle))
            etap = (n - 1j * k) * (n - 1j * k) / etas
            deltas = 2 * np.pi / lamda * self.entry_d_list[i] * etas
            deltap = 2 * np.pi / lamda * self.entry_d_list[i] * etas
            self.etas_list.append(etas)
            self.etap_list.append(etap)
            self.deltas_list.append(deltas)
            self.deltap_list.append(deltap)

        if self.subtype == 1:
            self.nsub = np.sqrt(11.67316 + 1 / lamda / lamda + 0.004482633 / (lamda * lamda - 1.108205 * 1.108205))
        elif self.subtype == 2:
            self.nsub = np.sqrt(5.68 + 1.53 * lamda * lamda / (lamda * lamda - 0.366))
        self.etasubs = np.sqrt(self.nsub * self.nsub - np.sin(self.angle) * np.sin(self.angle))
        self.etasubp = self.nsub * self.nsub / self.etasubs

        for i in range(0, len(self.layertype_list)):
            matrixs = np.matrix([[np.cos(self.deltas_list[len(self.layertype_list) - i - 1]), 1j * np.sin(self.deltas_list[len(self.layertype_list) - i - 1]) / self.etas_list[len(self.layertype_list) - i - 1]],
                                       [1j * self.etas_list[len(self.layertype_list) - i - 1] * np.sin(self.deltas_list[len(self.layertype_list) - i - 1]), np.cos(self.deltas_list[len(self.layertype_list) - i - 1])]])
            matrixp = np.matrix([[np.cos(self.deltap_list[len(self.layertype_list) - i - 1]), 1j * np.sin(self.deltap_list[len(self.layertype_list) - i - 1]) / self.etap_list[len(self.layertype_list) - i - 1]],
                                       [1j * self.etap_list[len(self.layertype_list) - i - 1] * np.sin(self.deltap_list[len(self.layertype_list) - i - 1]), np.cos(self.deltap_list[len(self.layertype_list) - i - 1])]])

            self.matrixs_list.append(matrixs)
            self.matrixp_list.append(matrixp)

        submatrixs = np.array([[1], [self.etasubs]])
        submatrixs.reshape(2, 1)
        submatrixp = np.array([[1], [self.etasubp]])
        submatrixp.reshape(2, 1)

        products = self.matrixs_list[0]

        for i in range(1, len(self.matrixs_list)):
            products = np.dot(products, self.matrixs_list[i])

        products = np.dot(products, submatrixs)
        Bs = products.item(0)
        Cs = products.item(1)

        productp = self.matrixp_list[0]

        for i in range(1, len(self.matrixp_list)):
            productp = np.dot(productp, self.matrixp_list[i])

        productp = np.dot(productp, submatrixp)

        Bp = productp.item(0)
        Cp = productp.item(1)

        Zs = self.eta0s * Bs + Cs
        Zp = self.eta0p * Bp + Cp
        Z2s = self.eta0s * Bs - Cs
        Z2p = self.eta0p * Bp - Cp

        Ztops = Bs * (Cs.conjugate()) - self.etasubs
        Ztopp = Bp * (Cp.conjugate()) - self.etasubp

        Ts = 4 * self.eta0s * self.etasubs / Zs / Zs.conjugate()
        Tp = 4 * self.eta0p * self.etasubp / Zp / Zp.conjugate()

        Rs = Z2s / Zs * ((Z2s / Zs).conjugate())
        Rp = Z2p / Zp * ((Z2p / Zp).conjugate())

        As = 4 * self.eta0s * Ztops.real / Zs / Zs.conjugate()
        Ap = 4 * self.eta0p * Ztopp.real / Zp / Zp.conjugate()

        transmission = (Ts + Tp) / 2 * 100 * self.scalefactor
        reflection = (Rs + Rp) / 2 * 100 * 1 + (Ts + Tp) / 2 * 100 * (1 - self.scalefactor)
        absorption = (As + Ap) / 2 * 100 * 1

        self.allresult.append(float(transmission))
        self.allresult.append(float(reflection))
        self.allresult.append(float(absorption))

        return self.allresult

    def cal_initialpara(self, x):
        self.T = 300
        self.A1 = 13.173 - 9.852 * x + 2.909 * x * x + 0.0001 * (300 - self.T)
        self.B1 = 0.83 - 0.246 * x - 0.0961 * x * x + 8 * 0.00001 * (300 - self.T)
        self.C1 = 6.706 - 14.437 * x + 8.531 * x * x + 7 * 0.00001 * (300 - self.T)
        self.D1 = 1.953 * 0.00001 - 0.00128 * x + 1.853 * 0.00001 * x * x
        self.T0 = 81.9
        self.E0 = -0.3424 + 1.838 * x
        self.sigma = 3.267 * np.power(10, 4) * (1 + x)
        self.alpha0 = np.exp(-18.88 + 3.61 * x)

    def load_n_file(self):
        self.wl_n_ZnSe = []
        self.n_ZnSe = []
        self.wl_k_ZnSe = []
        self.k_ZnSe = []
        self.wl_n_BaF2 = []
        self.n_BaF2 = []
        self.wl_k_BaF2 = []
        self.k_BaF2 = []
        self.wl_n_Ge = []
        self.n_Ge = []
        self.wl_k_Ge = []
        self.k_Ge = []
        self.wl_n_ZnS = []
        self.n_ZnS = []
        self.wl_k_ZnS = []
        self.k_ZnS = []

        os.chdir('/Users/apple/Dropbox/6.python/Projects/FTIR_fittingtool_v2/Refractive_Index')
        with open('ZnSe_n.csv', 'r') as f1:
            reader = csv.reader(f1, delimiter=',')
            for row in reader:
                try:
                    self.wl_n_ZnSe.append(float(row[0]))
                    self.n_ZnSe.append(float(row[1]))
                except ValueError:
                    pass
        with open('ZnSe_k.csv', 'r') as f2:
            reader = csv.reader(f2, delimiter=',')
            for row in reader:
                try:
                    self.wl_k_ZnSe.append(float(row[0]))
                    self.k_ZnSe.append(float(row[1]))
                except ValueError:
                    pass
        with open('BaF2_n.csv', 'r') as f3:
            reader = csv.reader(f3, delimiter=',')
            for row in reader:
                try:
                    self.wl_n_BaF2.append(float(row[0]))
                    self.n_BaF2.append(float(row[1]))
                except ValueError:
                    pass
        with open('BaF2_k.csv', 'r') as f4:
            reader = csv.reader(f4, delimiter=',')
            for row in reader:
                try:
                    self.wl_k_BaF2.append(float(row[0]))
                    self.k_BaF2.append(float(row[1]))
                except ValueError:
                    pass

        with open('Ge_n_293K.csv', 'r') as f5:
            reader = csv.reader(f5, delimiter=',')
            for row in reader:
                try:
                    self.wl_n_Ge.append(float(row[0]))
                    self.n_Ge.append(float(row[1]))
                except ValueError:
                    pass

        with open('ZnS_n.csv', 'r') as f7:
            reader = csv.reader(f7, delimiter=',')
            for row in reader:
                try:
                    self.wl_n_ZnS.append(float(row[0]))
                    self.n_ZnS.append(float(row[1]))
                except ValueError:
                    pass

        with open('ZnS_k.csv', 'r') as f8:
            reader = csv.reader(f8, delimiter=',')
            for row in reader:
                try:
                    self.wl_k_ZnS.append(float(row[0]))
                    self.k_ZnS.append(float(row[1]))
                except ValueError:
                    pass

    def cal_n(self, lamda, material):
        if material == "CdTe" or material == "MCT" or material == "SL":
            if lamda < 1.4 * self.C1:
                lamda = 1.4 * self.C1

            n = np.sqrt(self.A1 + self.B1 / (1 - (self.C1 / lamda) * (self.C1 / lamda)) + self.D1 * lamda * lamda)
            return n
        elif material == "ZnSe":
            for i in range(0, len(self.wl_n_ZnSe)):
                if self.wl_n_ZnSe[i + 1] > lamda >= self.wl_n_ZnSe[i]:
                    n = self.n_ZnSe[i]
                    return n

        elif material == "BaF2":
            for i in range(0, len(self.wl_n_BaF2)):
                if self.wl_n_BaF2[i + 1] > lamda >= self.wl_n_BaF2[i]:
                    n = self.n_BaF2[i]
                    return n

        elif material == "Ge":
            try:
                for i in range(0, len(self.wl_n_Ge)):
                    if self.wl_n_Ge[i + 1] > lamda >= self.wl_n_Ge[i]:
                        n = self.n_Ge[i]
                        return n
            except IndexError:
                return 4  # The Ge refractive index file does not have enough range.

        elif material == "ZnS":
            for i in range(0, len(self.wl_n_ZnS)):
                if self.wl_n_ZnS[i + 1] > lamda >= self.wl_n_ZnS[i]:
                    n = self.n_ZnS[i]
                    return n

        elif material == "Si":
            n = np.sqrt(11.67316 + 1 / lamda / lamda + 0.004482633 / (lamda * lamda - 1.108205 * 1.108205))
            return n
        elif material == "Air":
            n = 1
            return n

    def cal_k(self, lamda, material):
        k = 0
        if material == "CdTe" or material == "MCT" or material == "SL":
            return 0

        elif material == "ZnSe":
            for i in range(0, len(self.wl_k_ZnSe)):
                if self.wl_k_ZnSe[i + 1] > lamda >= self.wl_k_ZnSe[i]:
                    k = self.k_ZnSe[i]
                    return k

        elif material == "BaF2":
            for i in range(0, len(self.wl_k_BaF2)):
                if self.wl_k_BaF2[i + 1] > lamda >= self.wl_k_BaF2[i]:
                    k = self.k_BaF2[i]
                    return k

        elif material == "Ge":
            return 0

        elif material == "ZnS":
            for i in range(0, len(self.wl_k_ZnS)):
                if self.wl_k_ZnS[i + 1] > lamda >= self.wl_k_ZnS[i]:
                    k = self.k_ZnS[i]
                    return k

        elif material == "Si":
            return 0

        elif material == "Air":
            return 0

    def cal_R(self, x1, x2, lamda):
        if x1 != -1:
            self.cal_initialpara(x1)
            n1 = self.cal_n(lamda)
        else:
            n1 = 1
        if x2 != -1:
            self.cal_initialpara(x2)
            n2 = self.cal_n(lamda)
        else:
            n2 = 1
        self.R = ((n1 - n2) / (n1 + n2)) * ((n1 - n2) / (n1 + n2))
        return self.R

    def cal_n_array(self, lamda):
        n = np.sqrt(self.A1 + self.B1 / (1 - (self.C1 / lamda) * (self.C1 / lamda)) + self.D1 * lamda * lamda)
        return n

    def cal_alpha(self):
        self.alpha = self.alpha0 * math.exp(self.sigma * (self.E - self.E0) / (self.T + self.T0))
        return self.alpha

    def cal_alpha_all(self):
        z = 0
        xs = self.x
        s = self.sfactor
        deltaz = self.bufferd
        self.a11 = 0
        while z < self.layerd:
            x_cal = (1 - (xs + s * self.layerd)) / (1 + 4 * (z / deltaz) * (z / deltaz)) + xs + s * self.layerd - s * z
            self.cal_initialpara(x_cal)
            self.alphaMCT += self.cal_alpha()
            self.a11 += self.alphaMCT * 0.05
            z += 0.05
        self.a1 = math.exp(-self.a11)
        self.cal_initialpara(1)
        self.alphaCT = self.cal_alpha()
        self.a2 = math.exp(-self.alphaCT * self.subd)

    def cal_trans(self):
        self.R23 = self.R2 + self.R3 * (1 - self.R2) * (1 - self.R2) * self.a2 * self.a2 / (
        1 - self.R2 * self.R3 * self.a2 * self.a2)
        self.T23 = (1 - self.R2) * (1 - self.R3) * self.a2 / (1 - self.R2 * self.R3 * self.a2 * self.a2)

        self.T13 = (1 - self.R1) * (1 - self.H) * self.T23 * self.a1 \
                   / (1 - self.R1 * (1 - self.H) * self.R23 * self.a1 * self.a1)
        self.addlog(self.T13)
        return self.T13

    def cal_fringes_old(self, lamda):
        self.peakvalue = []
        self.n_list = []
        self.k_list = []
        self.etas_list = []
        self.etap_list = []
        self.deltas_list = []
        self.deltap_list = []
        self.matrixs_list = []
        self.matrixp_list = []

        self.eta0s = np.cos(self.angle)
        self.eta0p = 1 / np.cos(self.angle)

        for i in range(0, len(self.layertype_list)):
            if self.layertype_list[i] == "CdTe":
                self.cal_initialpara(1)
            else:
                self.cal_initialpara(self.entry_x_list[i])
            n = self.cal_n_array(lamda)
            k = 0
            self.n_list.append(n)
            self.k_list.append(k)
            etas = np.sqrt((n - 1j * k) * (n - 1j * k) - np.sin(self.angle) * np.sin(self.angle))
            etap = (n - 1j * k) * (n - 1j * k) / etas
            deltas = 2 * np.pi / lamda * self.entry_d_list[i] * etas
            deltap = 2 * np.pi / lamda * self.entry_d_list[i] * etas
            self.etas_list.append(etas)
            self.etap_list.append(etap)
            self.deltas_list.append(deltas)
            self.deltap_list.append(deltap)

        if self.subtype == 1:
            self.nsub = np.sqrt(11.67316 + 1 / lamda / lamda + 0.004482633 / (lamda * lamda - 1.108205 * 1.108205))
        elif self.subtype == 2:
            self.nsub = np.sqrt(5.68 + 1.53 * lamda * lamda / (lamda * lamda - 0.366))
        self.etasubs = np.sqrt(self.nsub * self.nsub - np.sin(self.angle) * np.sin(self.angle))
        self.etasubp = self.nsub * self.nsub / self.etasubs

        for l in range(0, len(lamda)):
            for i in range(0, len(self.layertype_list)):
                matrixs = np.matrix([[np.cos(self.deltas_list[len(self.layertype_list) - i - 1][l]), 1j * np.sin(self.deltas_list[len(self.layertype_list) - i - 1][l]) / self.etas_list[len(self.layertype_list) - i - 1][l]],
                                           [1j * self.etas_list[len(self.layertype_list) - i - 1][l] * np.sin(self.deltas_list[len(self.layertype_list) - i - 1][l]), np.cos(self.deltas_list[len(self.layertype_list) - i - 1][l])]])
                matrixp = np.matrix([[np.cos(self.deltap_list[len(self.layertype_list) - i - 1][l]), 1j * np.sin(self.deltap_list[len(self.layertype_list) - i - 1][l]) / self.etap_list[len(self.layertype_list) - i - 1][l]],
                                           [1j * self.etap_list[len(self.layertype_list) - i - 1][l] * np.sin(self.deltap_list[len(self.layertype_list) - i - 1][l]), np.cos(self.deltap_list[len(self.layertype_list) - i - 1][l])]])

                self.matrixs_list.append(matrixs)
                self.matrixp_list.append(matrixp)

            submatrixs = np.array([[1], [self.etasubs[l]]])
            submatrixs.reshape(2, 1)
            submatrixp = np.array([[1], [self.etasubp[l]]])
            submatrixp.reshape(2, 1)

            products = self.matrixs_list[0]

            for i in range(1, len(self.matrixs_list)):
                products = np.dot(products, self.matrixs_list[i])

            products = np.dot(products, submatrixs)
            Bs = products.item(0)
            Cs = products.item(1)

            productp = self.matrixp_list[0]

            for i in range(1, len(self.matrixp_list)):
                productp = np.dot(productp, self.matrixp_list[i])

            productp = np.dot(productp, submatrixp)

            Bp = productp.item(0)
            Cp = productp.item(1)

            Zs = self.eta0s * Bs + Cs
            Zp = self.eta0p * Bp + Cp

            Ts = 4 * self.eta0s * self.etasubs[l] / Zs / Zs.conjugate()
            Tp = 4 * self.eta0p * self.etasubp[l] / Zp / Zp.conjugate()

            self.peakvalue.append(float((Ts + Tp)) / 2 * 100 * self.scalefactor)

        print(self.peakvalue)

        return self.peakvalue

    def cal_absorption(self):
        basek = 0
        angle = self.angle
        self.absorptions = []
        numbercount = 0

        scalefactor = self.scalefactor
        based = self.based

        eta0s = np.cos(angle)
        eta0p = 1 / np.cos(angle)

        for i in range(0, len(self.wns)):
            lamda = 10000 / self.wns[i]
            trans = self.trans[i]

            numbercount += 1
            if numbercount == 1:
                percentage = (self.wns[i] - self.wns[0]) / (self.wns[len(self.wns) - 1] - self.wns[0]) * 100
                self.progress_var.set(percentage)
                self.wn_beingcalculated.set(self.wns[i])
                # timemin = 16 / 5810 * len(self.wns) * (1 - percentage / 100)
                # second1, minleft = math.modf(timemin)
                # secleft = second1 * 60
                # self.addlog('{0:.2f}% Complete. Estimated time remaining: {1:02d}:{2:02d} '.format(percentage, int(minleft), int(secleft)))
                numbercount = 0

            delta = 10
            fitsuccess = 0

            self.cal_initialpara(1)
            n1 = self.cal_n(lamda)
            k1 = 0
            eta1s = np.sqrt((n1 - 1j * k1) * (n1 - 1j * k1) - np.sin(angle) * np.sin(angle))
            eta1p = (n1 - 1j * k1) * (n1 - 1j * k1) / eta1s
            delta1s = 2 * np.pi / lamda * self.capd * eta1s
            delta1p = 2 * np.pi / lamda * self.capd * eta1s

            self.cal_initialpara(self.mctcapx)
            n2 = self.cal_n(lamda)
            k2 = 0
            eta2s = np.sqrt((n2 - 1j * k2) * (n2 - 1j * k2) - np.sin(angle) * np.sin(angle))
            eta2p = (n2 - 1j * k2) * (n2 - 1j * k2) / eta2s
            delta2s = 2 * np.pi / lamda * self.mctcapd * eta2s
            delta2p = 2 * np.pi / lamda * self.mctcapd * eta2s

            self.cal_initialpara(self.x)
            n3 = self.cal_n(lamda)

            self.cal_initialpara(self.bufferx)
            n4 = self.cal_n(lamda)
            k4 = 0
            eta4s = np.sqrt((n4 - 1j * k4) * (n4 - 1j * k4) - np.sin(angle) * np.sin(angle))
            eta4p = (n4 - 1j * k4) * (n4 - 1j * k4) / eta4s
            delta4s = 2 * np.pi / lamda * self.bufferd * eta4s
            delta4p = 2 * np.pi / lamda * self.bufferd * eta4s

            self.cal_initialpara(1)
            n5 = self.cal_n(lamda)
            k5 = 0
            eta5s = np.sqrt((n5 - 1j * k5) * (n5 - 1j * k5) - np.sin(angle) * np.sin(angle))
            eta5p = (n5 - 1j * k5) * (n5 - 1j * k5) / eta5s
            delta5s = 2 * np.pi / lamda * self.subd * eta5s
            delta5p = 2 * np.pi / lamda * self.subd * eta5s

            if self.subtype == 1:
                nsub = np.sqrt(11.67316 + 1 / lamda / lamda + 0.004482633 / (lamda * lamda - 1.108205 * 1.108205))
            elif self.subtype == 2:
                nsub = np.sqrt(5.68 + 1.53 * lamda * lamda / (lamda * lamda - 0.366))
            etasubs = np.sqrt(nsub * nsub - np.sin(angle) * np.sin(angle))
            etasubp = nsub * nsub / etasubs

            matrix1s = np.matrix([[np.cos(delta1s), 1j * np.sin(delta1s) / eta1s],
                                  [1j * eta1s * np.sin(delta1s), np.cos(delta1s)]])
            matrix1p = np.matrix([[np.cos(delta1p), 1j * np.sin(delta1p) / eta1p],
                                  [1j * eta1p * np.sin(delta1p), np.cos(delta1p)]])

            matrix2s = np.matrix([[np.cos(delta2s), 1j * np.sin(delta2s) / eta2s],
                                  [1j * eta2s * np.sin(delta2s), np.cos(delta2s)]])
            matrix2p = np.matrix([[np.cos(delta2p), 1j * np.sin(delta2p) / eta2p],
                                  [1j * eta2p * np.sin(delta2p), np.cos(delta2p)]])

            submatrixs = np.array([[1], [etasubs]])
            submatrixs.reshape(2, 1)
            submatrixp = np.array([[1], [etasubp]])
            submatrixp.reshape(2, 1)

            product1s = np.dot(matrix1s, matrix2s)
            product1p = np.dot(matrix1p, matrix2p)

            for k in np.arange(0, 1, 0.001):
                k3 = k
                eta3s = np.sqrt((n3 - 1j * k3) * (n3 - 1j * k3) - np.sin(angle) * np.sin(angle))
                eta3p = (n3 - 1j * k3) * (n3 - 1j * k3) / eta3s
                delta3s = 2 * np.pi / lamda * based * eta3s
                delta3p = 2 * np.pi / lamda * based * eta3s

                matrix3s = np.matrix([[np.cos(delta3s), 1j * np.sin(delta3s) / eta3s],
                                      [1j * eta3s * np.sin(delta3s), np.cos(delta3s)]])
                matrix3p = np.matrix([[np.cos(delta3p), 1j * np.sin(delta3p) / eta3p],
                                      [1j * eta3p * np.sin(delta3p), np.cos(delta3p)]])

                matrix4s = np.matrix([[np.cos(delta4s), 1j * np.sin(delta4s) / eta4s],
                                      [1j * eta4s * np.sin(delta4s), np.cos(delta4s)]])
                matrix4p = np.matrix([[np.cos(delta4p), 1j * np.sin(delta4p) / eta4p],
                                      [1j * eta4p * np.sin(delta4p), np.cos(delta4p)]])

                matrix5s = np.matrix([[np.cos(delta5s), 1j * np.sin(delta5s) / eta5s],
                                      [1j * eta5s * np.sin(delta5s), np.cos(delta5s)]])
                matrix5p = np.matrix([[np.cos(delta5p), 1j * np.sin(delta5p) / eta5p],
                                      [1j * eta5p * np.sin(delta5p), np.cos(delta5p)]])

                product2s = np.dot(product1s, matrix3s)
                product3s = np.dot(product2s, matrix4s)
                product4s = np.dot(product3s, matrix5s)
                product5s = np.dot(product4s, submatrixs)

                Bs = product5s.item(0)
                Cs = product5s.item(1)

                product2p = np.dot(product1p, matrix3p)
                product3p = np.dot(product2p, matrix4p)
                product4p = np.dot(product3p, matrix5p)
                product5p = np.dot(product4p, submatrixp)
                Bp = product5p.item(0)
                Cp = product5p.item(1)

                Zs = eta0s * Bs + Cs
                Zp = eta0p * Bp + Cp

                Ts = 4 * eta0s * etasubs / Zs / Zs.conjugate()
                Tp = 4 * eta0p * etasubp / Zp / Zp.conjugate()

                peakvalue = (Ts + Tp) / 2 * 100 * scalefactor

                if abs(peakvalue - trans) <= delta:
                    fitsuccess = 1
                    delta = abs(peakvalue - trans)
                    basek = k

            if fitsuccess == 1:
                ab = 4 * np.pi * basek / lamda * 10000
                self.absorptions.append(ab)
            else:
                self.addlog('Fitting failed at wavenumber = {}cm-1'.format(self.wns[i]))
                self.absorptions.append(0)

        self.addlog('Fitting complete!')
        return self.absorptions

    def cal_absorption_single(self, wn):
        basek = 0
        i = 0
        angle = self.angle

        scalefactor = self.scalefactor
        based = self.based

        eta0s = np.cos(angle)
        eta0p = 1 / np.cos(angle)

        while True:
            if self.wns[i + 1] > wn >= self.wns[i]:
                break
            else:
                i += 1

        lamda = 10000 / self.wns[i]
        trans = self.trans[i]

        delta = 10
        fitsuccess = 0

        self.cal_initialpara(1)
        n1 = self.cal_n(lamda)
        k1 = 0
        eta1s = np.sqrt((n1 - 1j * k1) * (n1 - 1j * k1) - np.sin(angle) * np.sin(angle))
        eta1p = (n1 - 1j * k1) * (n1 - 1j * k1) / eta1s
        delta1s = 2 * np.pi / lamda * self.capd * eta1s
        delta1p = 2 * np.pi / lamda * self.capd * eta1s

        self.cal_initialpara(self.mctcapx)
        n2 = self.cal_n(lamda)
        k2 = 0
        eta2s = np.sqrt((n2 - 1j * k2) * (n2 - 1j * k2) - np.sin(angle) * np.sin(angle))
        eta2p = (n2 - 1j * k2) * (n2 - 1j * k2) / eta2s
        delta2s = 2 * np.pi / lamda * self.mctcapd * eta2s
        delta2p = 2 * np.pi / lamda * self.mctcapd * eta2s

        self.cal_initialpara(self.x)
        n3 = self.cal_n(lamda)

        self.cal_initialpara(self.bufferx)
        n4 = self.cal_n(lamda)
        k4 = 0
        eta4s = np.sqrt((n4 - 1j * k4) * (n4 - 1j * k4) - np.sin(angle) * np.sin(angle))
        eta4p = (n4 - 1j * k4) * (n4 - 1j * k4) / eta4s
        delta4s = 2 * np.pi / lamda * self.bufferd * eta4s
        delta4p = 2 * np.pi / lamda * self.bufferd * eta4s

        self.cal_initialpara(1)
        n5 = self.cal_n(lamda)
        k5 = 0
        eta5s = np.sqrt((n5 - 1j * k5) * (n5 - 1j * k5) - np.sin(angle) * np.sin(angle))
        eta5p = (n5 - 1j * k5) * (n5 - 1j * k5) / eta5s
        delta5s = 2 * np.pi / lamda * self.subd * eta5s
        delta5p = 2 * np.pi / lamda * self.subd * eta5s

        if self.subtype == 1:
            nsub = np.sqrt(11.67316 + 1 / lamda / lamda + 0.004482633 / (lamda * lamda - 1.108205 * 1.108205))
        elif self.subtype == 2:
            nsub = np.sqrt(5.68 + 1.53 * lamda * lamda / (lamda * lamda - 0.366))
        etasubs = np.sqrt(nsub * nsub - np.sin(angle) * np.sin(angle))
        etasubp = nsub * nsub / etasubs

        matrix1s = np.matrix([[np.cos(delta1s), 1j * np.sin(delta1s) / eta1s],
                              [1j * eta1s * np.sin(delta1s), np.cos(delta1s)]])
        matrix1p = np.matrix([[np.cos(delta1p), 1j * np.sin(delta1p) / eta1p],
                              [1j * eta1p * np.sin(delta1p), np.cos(delta1p)]])

        matrix2s = np.matrix([[np.cos(delta2s), 1j * np.sin(delta2s) / eta2s],
                              [1j * eta2s * np.sin(delta2s), np.cos(delta2s)]])
        matrix2p = np.matrix([[np.cos(delta2p), 1j * np.sin(delta2p) / eta2p],
                              [1j * eta2p * np.sin(delta2p), np.cos(delta2p)]])

        submatrixs = np.array([[1], [etasubs]])
        submatrixs.reshape(2, 1)
        submatrixp = np.array([[1], [etasubp]])
        submatrixp.reshape(2, 1)

        product1s = np.dot(matrix1s, matrix2s)
        product1p = np.dot(matrix1p, matrix2p)

        for k in np.arange(0, 1, 0.001):
            k3 = k
            eta3s = np.sqrt((n3 - 1j * k3) * (n3 - 1j * k3) - np.sin(angle) * np.sin(angle))
            eta3p = (n3 - 1j * k3) * (n3 - 1j * k3) / eta3s
            delta3s = 2 * np.pi / lamda * based * eta3s
            delta3p = 2 * np.pi / lamda * based * eta3s

            matrix3s = np.matrix([[np.cos(delta3s), 1j * np.sin(delta3s) / eta3s],
                                  [1j * eta3s * np.sin(delta3s), np.cos(delta3s)]])
            matrix3p = np.matrix([[np.cos(delta3p), 1j * np.sin(delta3p) / eta3p],
                                  [1j * eta3p * np.sin(delta3p), np.cos(delta3p)]])

            matrix4s = np.matrix([[np.cos(delta4s), 1j * np.sin(delta4s) / eta4s],
                                  [1j * eta4s * np.sin(delta4s), np.cos(delta4s)]])
            matrix4p = np.matrix([[np.cos(delta4p), 1j * np.sin(delta4p) / eta4p],
                                  [1j * eta4p * np.sin(delta4p), np.cos(delta4p)]])

            matrix5s = np.matrix([[np.cos(delta5s), 1j * np.sin(delta5s) / eta5s],
                                  [1j * eta5s * np.sin(delta5s), np.cos(delta5s)]])
            matrix5p = np.matrix([[np.cos(delta5p), 1j * np.sin(delta5p) / eta5p],
                                  [1j * eta5p * np.sin(delta5p), np.cos(delta5p)]])

            product2s = np.dot(product1s, matrix3s)
            product3s = np.dot(product2s, matrix4s)
            product4s = np.dot(product3s, matrix5s)
            product5s = np.dot(product4s, submatrixs)

            Bs = product5s.item(0)
            Cs = product5s.item(1)

            product2p = np.dot(product1p, matrix3p)
            product3p = np.dot(product2p, matrix4p)
            product4p = np.dot(product3p, matrix5p)
            product5p = np.dot(product4p, submatrixp)
            Bp = product5p.item(0)
            Cp = product5p.item(1)

            Zs = eta0s * Bs + Cs
            Zp = eta0p * Bp + Cp

            Ts = 4 * eta0s * etasubs / Zs / Zs.conjugate()
            Tp = 4 * eta0p * etasubp / Zp / Zp.conjugate()

            peakvalue = (Ts + Tp) / 2 * 100 * scalefactor

            if abs(peakvalue - trans) <= delta:
                fitsuccess = 1
                delta = abs(peakvalue - trans)
                basek = k

        if fitsuccess == 1:
            ab = 4 * np.pi * basek / lamda * 10000
            return ab
        else:
            return 0

    def returntrans(self):
        return self.transs

    def returnpeakvalues(self):
        return self.peakvalues

    def returnreflections(self):
        return self.reflections

    def returnabsorptions(self):
        return self.absorptions

    def addlog(self, string):
        self.listbox.insert(END, string)
        self.listbox.yview(END)


class cal_MCT_a:
    def __init__(self, x, wavenumbers, fittype):
        self.x = float(x)
        self.wavenumbers = wavenumbers
        self.absorptions = []
        self.T = 300
        self.fittype = fittype

        if self.fittype == "Chu":
            self.a0 = np.exp(-18.5 + 45.68 * self.x)
            self.E0 = -0.355 + 1.77 * self.x
            self.ag = -65 + 1.88 * self.T + (8694 - 10.31 * self.T) * self.x
            self.Eg = -0.295 + 1.87 * self.x - 0.28 * self.x * self.x + \
                      (6 - 14 * self.x + 3 * self.x * self.x) * 0.0001 * self.T + 0.35 * self.x * self.x * self.x * self.x
            self.delta_kT = (np.log(self.ag) - np.log(self.a0)) / (self.Eg - self.E0)
            self.beta = -1 + 0.083 * self.T + (21 - 0.13 * self.T) * self.x
        elif self.fittype == "Schacham and Finkman":
            self.a0 = np.exp(-18.88 + 53.61 * self.x)
            self.E0 = -0.3424 + 1.838 * self.x + 0.148 * self.x * self.x * self.x * self.x
            self.Eg = self.E0 + (0.0629 + 0.000768 * self.T) * (1 - 2.14 * self.x) / (1 + self.x)
            self.delta_kT = 32670 * (1 + self.x) / (self.T + 81.9)
            self.beta = 210900 * np.sqrt((1 + self.x) / (81.9 + self.T))
        elif self.fittype == "Yong":
            self.W = 0.013 # not sure
            self.A = 1 / np.power(10, 17) # not sure
            self.P = 8 / np.power(10, 8)
            self.s = np.sqrt(2 * self.P * self.P / 3)
            self.B = self.A / np.power(np.pi, 2) / np.power(self.s, 3)
            self.E0 = -0.295 + 1.87 * self.x - 0.28 * self.x * self.x + \
                      (6 - 14 * self.x + 3 * self.x * self.x) * 0.0001 * self.T + 0.35 * self.x * self.x * self.x * self.x

            print(self.E0)
            # Here E0 is the same equation as Eg in Chu's fomula.
            self.Eg = self.E0 - self.W / 2
            self.b = self.Eg / 2

        self.ab = 0

        self.cal_all()

    def cal_Urbach(self, energy):
        if self.fittype == "Chu" or self.fittype == "Schacham and Finkman":
            self.ab = self.a0 * np.exp(self.delta_kT * (energy - self.E0))
        elif self.fittype == "Yong":
            self.ab = self.B / self.E0 * \
                      ((self.W / 2 + self.b) * np.sqrt((self.W / 2 + self.b) * (self.W / 2 + self.b) - self.b * self.b)
                       + 1 / 8 * (self.W / 2 + 2 * self.b) * np.sqrt((self.W / 2 + 2 * self.b) * (self.W / 2 + 2 * self.b) - 4 * self.b * self.b)) * np.exp(energy / self.W)

        return self.ab

    def cal_Kane(self, energy):
        if self.fittype == "Chu":
            self.ab = self.ag * np.exp(np.sqrt(self.beta * (energy - self.Eg)))
        elif self.fittype == "Schacham and Finkman":
            self.ab = self.beta * np.sqrt(energy - self.Eg)
        elif self.fittype == "Yong":
            self.ab = self.B / energy * ((energy - self.Eg + self.b) * np.sqrt((energy - self.Eg + self.b) * (energy - self.Eg + self.b) - self.b * self.b)
                                         + 1 / 8 * (energy - self.Eg + 2 * self.b) * np.sqrt((energy - self.Eg + 2 * self.b) * (energy - self.Eg + 2 * self.b) - 4 * self.b * self.b))
        return self.ab

    def cal_all(self):
        self.absorptions = []
        for wavenumber in self.wavenumbers:
            wl = 10000 / wavenumber
            E = 4.13566743 * 3 / 10 / wl      # Here E is in unit of eV.
            if self.fittype == "Chu" or self.fittype == "Schacham and Finkman":
                if E <= self.Eg:
                    self.absorptions.append(self.cal_Urbach(E))
                else:
                    self.absorptions.append(self.cal_Kane(E))
            if self.fittype == "Yong":
                if E <= self.Eg + self.W/2:
                    self.absorptions.append(self.cal_Urbach(E))
                else:
                    self.absorptions.append(self.cal_Kane(E))

    def return_absorptions(self):
        return self.absorptions


class ThreadedTask1(threading.Thread):
    def __init__(self, queue_1, composition, layerd, wavenumbers, transmissions, capd, mctcapd, based, bufferd, subd,
                 scalefactor, sfactor, factorh, angle, mctcapx, bufferx, maxtrans, subtype, fittype, listbox,
                 progress_var, wn_beingcalculated):
        threading.Thread.__init__(self)
        self.queue = queue_1
        self.composition = composition
        self.layerd = layerd
        self.wavenumbers = wavenumbers
        self.transmissions = transmissions
        self.capd = capd
        self.mctcapd = mctcapd
        self.based = based
        self.bufferd = bufferd
        self.subd = subd
        self.scalefactor = scalefactor
        self.sfactor = sfactor
        self.factorh = factorh
        self.angle = angle
        self.mctcapx = mctcapx
        self.bufferx = bufferx
        self.maxtrans = maxtrans
        self.subtype = subtype
        self.fittype = fittype
        self.listbox = listbox
        self.progress_var = progress_var
        self.wn_beingcalculated = wn_beingcalculated

    def run(self):
        fitobject = FIT_FTIR(self.composition, self.layerd, self.wavenumbers, self.transmissions, self.capd,
                             self.mctcapd,
                             self.based, self.bufferd, self.subd, self.scalefactor, self.sfactor, self.factorh,
                             self.angle,
                             self.mctcapx, self.bufferx, self.maxtrans, self.subtype, self.fittype, self.listbox,
                             self.progress_var, self.wn_beingcalculated)
        self.queue.put(fitobject.cal_absorption())


class ThreadedTask_fringes(threading.Thread):
    def __init__(self, queue_1, inital_CdTe, inital_HgTe, entry_d_list_initial, layernumber, wavenumbers_cut, trans_cut,
                 subd, layertype_list, entry_x_list, entry_d_list, checklayer_list, scalefactor, angle,
                 subtype, fittype, listbox, progress_var, wn_beingcalculated):
        threading.Thread.__init__(self)
        self.queue = queue_1
        self.inital_CdTe = inital_CdTe
        self.inital_HgTe = inital_HgTe
        self.entry_d_list_initial = entry_d_list_initial
        self.layernumber = layernumber
        self.wavenumbers_cut = wavenumbers_cut
        self.trans_cut = trans_cut
        self.scalefactor = scalefactor
        self.angle = angle
        self.subd = subd
        self.layertype_list, self.entry_x_list, self.entry_d_list, self.checklayer_list \
            = layertype_list, entry_x_list, entry_d_list, checklayer_list

        self.peakvalues_fit = []
        self.MSE = 0
        self.smallest_MSE = 100000000
        self.best_CdTe_offset = 0
        self.best_HgTe_offset = 0

        self.subtype = subtype
        self.fittype = fittype
        self.listbox = listbox
        self.progress_var = progress_var
        self.wn_beingcalculated = wn_beingcalculated

    def run(self):
        CdTe_fitrange = 10
        HgTe_fitrange = 5
        for CdTe_offset in np.arange(self.inital_CdTe - CdTe_fitrange, self.inital_CdTe + CdTe_fitrange, 1):
            for HgTe_offset in np.arange(self.inital_HgTe - HgTe_fitrange, self.inital_HgTe + HgTe_fitrange, 1):
                self.progress_var.set((CdTe_offset - self.inital_CdTe + CdTe_fitrange) / CdTe_fitrange / 2 * 100 +
                                      (HgTe_offset - self.inital_HgTe + HgTe_fitrange) / HgTe_fitrange / 2 / CdTe_fitrange / 2 * 100)
                for i in range(0, self.layernumber):
                    if int(self.checklayer_list[i]) == 1:
                        if self.layertype_list[i] == "CdTe":
                            new_d = float(self.entry_d_list_initial[i]) * (1 + 0.01 * CdTe_offset)
                            self.entry_d_list[i] = new_d
                        elif self.layertype_list[i] == "MCT" or self.layertype_list[i] == "SL":
                            new_d = float(self.entry_d_list_initial[i]) * float(self.entry_x_list[i]) * (1 + 0.01 * CdTe_offset) \
                                    + float(self.entry_d_list_initial[i]) * (1 - float(self.entry_x_list[i])) * (1 + 0.01 * HgTe_offset)
                            self.entry_d_list[i] = new_d

                fitobject = FIT_FTIR(self.wavenumbers_cut, self.trans_cut, self.subd, self.layertype_list, self.entry_x_list,
                                     self.entry_d_list, self.checklayer_list, self.scalefactor, self.angle,
                                     CdTe_offset, HgTe_offset, self.subtype, 2, self.listbox, self.progress_var, self.wn_beingcalculated)
                self.peakvalues_fit = fitobject.returnpeakvalues()

                self.MSE = 0

                for i in range(0, len(self.trans_cut)):
                    self.MSE += (self.peakvalues_fit[i] - self.trans_cut[i]) * \
                                (self.peakvalues_fit[i] - self.trans_cut[i]) / len(self.trans_cut)

                if self.MSE <= self.smallest_MSE:
                    self.smallest_MSE = self.MSE
                    self.best_CdTe_offset = CdTe_offset
                    self.best_HgTe_offset = HgTe_offset

        result = [self.best_CdTe_offset, self.best_HgTe_offset, self.smallest_MSE]

        self.queue.put(result)


class FTIR_fittingtool_GUI(Frame):
    def __init__(self, root, masterroot, listbox, statusbar, status1, status2):
        super().__init__(root, width=cross_platform_config.config.FRAME_WIDTH, bg='#2b2b2b')
        self.root = root
        self.masterroot = masterroot
        self.listbox = listbox
        self.statusbar = statusbar
        self.status1 = status1
        self.status2 = status2
        self.wavenumbers = []
        self.wavenumbers_cut = []
        self.transmissions = []
        self.trans_cut = []
        self.absorptions = []
        self.filename = ''
        self.numberofdata = 0
        self.numberofdata2 = 0
        self.colororders = ['blue', 'green', 'cyan', 'magenta', 'yellow', 'black']
        self.colororders2 = ['red', 'green', 'cyan', 'magenta', 'yellow', 'black', 'red', 'green', 'cyan', 'magenta',
                             'yellow', 'black']
        self.trans = 0
        self.wavenumber = 0
        self.wavelength = 0
        self.energy = 0
        self.composition = 0
        self.Temp = 300
        self.fittype = IntVar()
        self.lowercut = 400
        self.highercut = 6000
        self.transcut = 70
        self.transcutlow = 0
        self.transmissions_fit = []
        self.peakvalues_fit = []
        self.reflections_fit = []
        self.absorptions_fit = []
        self.fitline = None
        self.MSE = 0
        self.MSE_new = 0
        self.fittedthickness = 0
        self.subtype = 1
        self.progress_var = DoubleVar()
        self.text = ''
        self.text2 = ''
        self.wn_beingcalculated = DoubleVar()

        self.osdir = ''

        self.available_materials = ["CdTe", "MCT", "SL", "Si", "ZnSe", "BaF2", "Ge", "ZnS", "Air"]

        self.displayreflection, self.displayabsorption = 0, 0

        self.COLUMN0_WIDTH = 4
        self.COLUMN1_WIDTH = 3
        self.COLUMN2_WIDTH = 8

        self.frame0 = Frame(self, width=1000, height=20, bg='#262626', bd=1)
        self.frame0.pack(side=TOP, fill=X, expand=True)
        self.frame0.pack_propagate(0)

        buttonsettings = Button(self.frame0, text="Settings",
                            command=self.settings, highlightbackground='#262626', width=7)
        buttonsettings.pack(side=LEFT)

        buttonopen = Button(self.frame0, text="Open(⌘+O)",
                            command=self.openfromfile, highlightbackground='#262626', width=9)
        buttonopen.pack(side=LEFT)

        buttonsave= Button(self.frame0, text="Save result",
                            command=self.savetofile, highlightbackground='#262626', width=9)
        buttonsave.pack(side=LEFT)

        buttonclear = Button(self.frame0, text="Clear(⌘+C)",
                             command=self.clearalldata, highlightbackground='#262626', width=9)
        buttonclear.pack(side=LEFT)

        buttonload2 = Button(self.frame0, text="Load Structure(⌘+L)",
                           command=self.load_structure, highlightbackground='#262626', width=15)
        buttonload2.pack(side=LEFT)

        buttonsave2 = Button(self.frame0, text="Save Structure",
                           command=self.save_structure, highlightbackground='#262626', width=12)
        buttonsave2.pack(side=LEFT)

        buttonsave = Button(self.frame0, text="Show Fringes (⌘+S)",
                            command=self.show_fringes, highlightbackground='#262626', width=15)
        buttonsave.pack(side=RIGHT)

        buttonfringes = Button(self.frame0, text="FIT Fringes(⌘+G)",
                               command=self.fit_fringes, highlightbackground='#262626', width=13)
        buttonfringes.pack(side=RIGHT)

        buttoncal = Button(self.frame0, text="Cal a(⌘+A)",
                           command=self.cal_absorption, highlightbackground='#262626', width=9)
        buttoncal.pack(side=RIGHT)

        buttonmct = Button(self.frame0, text="MCT a",
                           command=self.cal_MCT_absorption, highlightbackground='#262626', width=5)
        buttonmct.pack(side=RIGHT)

        self.filepath = Label(self.frame0, text="", bg='#262626', fg="#a9b7c6", width=23)
        self.filepath.pack(side=LEFT, fill=X)

        if _platform == "win32" or _platform == "win64":
            buttonopen.config(text='Open(Ct+O)')
            buttonclear.config(text='Clear(Ct+C)')
            buttoncal.config(text='Cal a(Ct+A)')
            buttonload2.config(text='Load Structure(Ct+L)')
            buttonfringes.config(text='FIT Fringes(Ct+G)')
            buttonsave.config(text='Show Fringes (Ct+S)')

        self.frame3 = Frame(self, width=150, bg='#2b2b2b')
        self.frame3.pack(side=RIGHT, fill=Y, expand=True)
        self.frame3.pack_propagate(0)

        LABEL_WIDTH = 13

        # seperateline = Label(self.frame3, text='-'*40, anchor=E, fg="#a9b7c6", bg='#2b2b2b')
        # seperateline.grid(row=9, column=0, columnspan=4)

        def change_sub(*args):
            if self.varsub.get() == "CdTe/Si":
                self.subtype = 1
            elif self.varsub.get() == "CdZnTe":
                self.subtype = 2

        Label(self.frame3, text='Layer:', fg="#a9b7c6", bg='#2b2b2b', width=self.COLUMN0_WIDTH).grid(row=10, column=0, sticky=E)
        Label(self.frame3, text='x:', fg="#a9b7c6", bg='#2b2b2b', width=self.COLUMN1_WIDTH + 1).grid(row=10, column=1, sticky=E)
        Label(self.frame3, text='d:', fg="#a9b7c6", bg='#2b2b2b', width=self.COLUMN2_WIDTH).grid(row=10, column=2, sticky=E)

        self.varsub = StringVar(self.frame3)
        self.varsub.set("Si")  # initial value
        self.varsub.trace("w", change_sub)
        suboption1 = OptionMenu(self.frame3, self.varsub, "Si", "CdZnTe")
        suboption1.config(bg="#2b2b2b", width=self.COLUMN0_WIDTH, anchor=E)
        suboption1.grid(row=11, column=0, sticky=W+E)

        self.entry_d_0 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_0.grid(row=11, column=2)
        self.entry_d_0.insert(0, '500')

        Label(self.frame3, text='  ', fg="#a9b7c6", bg='#2b2b2b', width=self.COLUMN0_WIDTH).grid(row=11, column=3, sticky=E)

        self.layernumber = 0
        self.layertypevar1, self.layertypevar2, self.layertypevar3, self.layertypevar4, self.layertypevar5, \
        self.layertypevar6, self.layertypevar7, self.layertypevar8, self.layertypevar9 , self.layertypevar10, \
        self.layertypevar11, self.layertypevar12, self.layertypevar13, self.layertypevar14, self.layertypevar15, \
        self.layertypevar16, self.layertypevar17, self.layertypevar18 \
            = StringVar(self.frame3), StringVar(self.frame3), StringVar(self.frame3), StringVar(self.frame3), \
              StringVar(self.frame3), StringVar(self.frame3), StringVar(self.frame3), StringVar(self.frame3), \
              StringVar(self.frame3), StringVar(self.frame3), StringVar(self.frame3), StringVar(self.frame3), \
              StringVar(self.frame3), StringVar(self.frame3), StringVar(self.frame3), StringVar(self.frame3), \
              StringVar(self.frame3), StringVar(self.frame3)
        self.checklayer1, self.checklayer2, self.checklayer3, self.checklayer4, self.checklayer5, \
        self.checklayer6, self.checklayer7, self.checklayer8, self.checklayer9, self.checklayer10, \
        self.checklayer11, self.checklayer12, self.checklayer13, self.checklayer14, self.checklayer15, \
        self.checklayer16, self.checklayer17, self.checklayer18 \
            = IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), \
              IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), IntVar()

        self.entry_x_1 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN1_WIDTH)
        self.entry_x_2 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN1_WIDTH)
        self.entry_x_3 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN1_WIDTH)
        self.entry_x_4 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN1_WIDTH)
        self.entry_x_5 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN1_WIDTH)
        self.entry_x_6 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN1_WIDTH)
        self.entry_x_7 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN1_WIDTH)
        self.entry_x_8 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN1_WIDTH)
        self.entry_x_9 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN1_WIDTH)
        self.entry_x_10 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN1_WIDTH)
        self.entry_x_11 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN1_WIDTH)
        self.entry_x_12 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN1_WIDTH)
        self.entry_x_13 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN1_WIDTH)
        self.entry_x_14 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN1_WIDTH)
        self.entry_x_15 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN1_WIDTH)
        self.entry_x_16 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN1_WIDTH)
        self.entry_x_17 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN1_WIDTH)
        self.entry_x_18 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN1_WIDTH)

        self.entry_d_1 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_2 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_3 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_4 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_5 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_6 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_7 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_8 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_9 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_10 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_11 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_12 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_13 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_14 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_15 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_16 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_17 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_18 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)

        self.layertype_list, self.entry_x_list, self.entry_d_list, self.checklayer_list = [], [], [], []

        def add_layer_on_top():
            self.layernumber += 1

            self.buttonaddlayer.grid_forget()
            if self.layernumber == 1:
                getattr(self, "layertypevar{}".format(self.layernumber)).set("CdTe")
            else:
                getattr(self, "layertypevar{}".format(self.layernumber)).set("MCT")
            getattr(self, "layertypevar{}".format(self.layernumber)).trace("w", change_sub)
            layertypeoption1 = OptionMenu(self.frame3, getattr(self, "layertypevar{}".format(self.layernumber)), *self.available_materials)
            layertypeoption1.config(bg="#2b2b2b", width=self.COLUMN0_WIDTH, anchor=E)
            layertypeoption1.grid(row=11 + self.layernumber, column=0, sticky=W+E)

            getattr(self, "entry_x_{}".format(self.layernumber)).grid(row=11 + self.layernumber, column=1, sticky=W+E)
            if self.layernumber == 1:
                getattr(self, "entry_x_{}".format(self.layernumber)).insert(0, '1.00')
            else:
                getattr(self, "entry_x_{}".format(self.layernumber)).insert(0, '0.30')

            getattr(self, "entry_d_{}".format(self.layernumber)).grid(row=11 + self.layernumber, column=2)
            if self.layernumber == 1:
                getattr(self, "entry_d_{}".format(self.layernumber)).insert(0, '11.0')
            else:
                getattr(self, "entry_d_{}".format(self.layernumber)).insert(0, '0.0')

            checkbox1 = Checkbutton(self.frame3, text="", variable=getattr(self, "checklayer{}".format(self.layernumber)), bg='#2b2b2b')
            checkbox1.grid(row=11 + self.layernumber, column=3, sticky=W)
            checkbox1.select()

            if self.layernumber < 16:
                self.buttonaddlayer.grid(row=12 + self.layernumber, column=0, columnspan=4, sticky=W+E)

        self.buttonaddlayer = Button(self.frame3, text="Add layer on top", command=add_layer_on_top, highlightbackground='#2b2b2b', width=12)
        self.buttonaddlayer.grid(row=12, column=0, columnspan=4, sticky=W+E)

        self.frame4 = Frame(self, width=850, bg='#2b2b2b')
        self.frame4.pack(side=TOP, fill=X, expand=True)
        self.frame4.pack_propagate(0)

        Label(self.frame4, text='WaveNum(cm-1):',
              fg="#a9b7c6", bg='#2b2b2b', width=LABEL_WIDTH, anchor=E).grid(row=0, column=0, columnspan=1, sticky=E)
        Label(self.frame4, text='Wavelength(um):',
              fg="#a9b7c6", bg='#2b2b2b', width=LABEL_WIDTH, anchor=E).grid(row=0, column=2, columnspan=1, sticky=E)
        Label(self.frame4, text='Energy(meV):',
              fg="#a9b7c6", bg='#2b2b2b', width=LABEL_WIDTH - 2, anchor=E).grid(row=0, column=4, columnspan=1, sticky=E)
        Label(self.frame4, text='MCT x:',
              fg="#a9b7c6", bg='#2b2b2b', width=LABEL_WIDTH - 4, anchor=E).grid(row=0, column=6, columnspan=1, sticky=E)

        self.wavenumber1 = Label(self.frame4, text='{}'.format(self.wavenumber), fg="#a9b7c6", bg='#2b2b2b', width=8,
                                 anchor=W)
        self.wavenumber1.grid(row=0, column=1, columnspan=1, sticky=E)
        self.wavelength1 = Label(self.frame4, text='{}'.format(self.wavelength), fg="#a9b7c6", bg='#2b2b2b', width=8,
                                 anchor=W)
        self.wavelength1.grid(row=0, column=3, columnspan=1, sticky=E)
        self.energy1 = Label(self.frame4, text='{}'.format(self.energy), fg="#a9b7c6", bg='#2b2b2b', width=8, anchor=W)
        self.energy1.grid(row=0, column=5, columnspan=1, sticky=E)
        self.composition1 = Label(self.frame4, text='{}'.format(self.composition), fg="#a9b7c6", bg='#2b2b2b', width=6,
                                  anchor=W)
        self.composition1.grid(row=0, column=7, columnspan=1, sticky=E)

        self.frame2 = Frame(self, width=850, bg='#2b2b2b')
        self.frame2.pack(side=BOTTOM, fill=X, expand=True)
        self.frame2.pack_propagate(0)

        label_31 = Label(self.frame2, text='X Lower Cut:', width=LABEL_WIDTH - 3, anchor=E, fg="#a9b7c6", bg='#2b2b2b')
        label_32 = Label(self.frame2, text='X Higher Cut:', width=LABEL_WIDTH - 3, anchor=E, fg="#a9b7c6", bg='#2b2b2b')
        label_33 = Label(self.frame2, text='Y Lower Cut:', width=LABEL_WIDTH - 3, anchor=E, fg="#a9b7c6", bg='#2b2b2b')
        label_332 = Label(self.frame2, text='Y Higher Cut:', width=LABEL_WIDTH - 3, anchor=E, fg="#a9b7c6", bg='#2b2b2b')
        self.entry_31 = Entry(self.frame2, highlightbackground='#2b2b2b', width=4)
        self.entry_32 = Entry(self.frame2, highlightbackground='#2b2b2b', width=4)
        self.entry_33 = Entry(self.frame2, highlightbackground='#2b2b2b', width=4)
        self.entry_332 = Entry(self.frame2, highlightbackground='#2b2b2b', width=4)
        label_31.grid(row=0, column=0, columnspan=1, sticky=E)
        label_32.grid(row=1, column=0, columnspan=1, sticky=E)
        label_33.grid(row=0, column=3, columnspan=1, sticky=E)
        label_332.grid(row=1, column=3, columnspan=1, sticky=E)
        self.entry_31.grid(row=0, column=1)
        self.entry_32.grid(row=1, column=1)
        self.entry_33.grid(row=0, column=4)
        self.entry_332.grid(row=1, column=4)
        self.entry_31.insert(0, self.lowercut)
        self.entry_32.insert(0, self.highercut)
        self.entry_33.insert(0, self.transcutlow)
        self.entry_332.insert(0, self.transcut)

        def getbutton31():
            self.entry_31.delete(0, END)
            self.entry_31.insert(0, "%.4f" % self.xclick)

        button31 = Button(self.frame2, text="Get",
                          command=getbutton31, highlightbackground='#2b2b2b', anchor=W, width=3)
        button31.grid(row=0, column=2, sticky=W)

        def getbutton32():
            self.entry_32.delete(0, END)
            self.entry_32.insert(0, "%.4f" % self.xclick)

        button32 = Button(self.frame2, text="Get",
                          command=getbutton32, highlightbackground='#2b2b2b', anchor=W, width=3)
        button32.grid(row=1, column=2, sticky=W)

        def getbutton33():
            self.entry_33.delete(0, END)
            self.entry_33.insert(0, "%.4f" % self.yclick)

        button33 = Button(self.frame2, text="Get",
                          command=getbutton33, highlightbackground='#2b2b2b', anchor=W, width=3)
        button33.grid(row=0, column=5, sticky=W)

        def getbutton332():
            self.entry_332.delete(0, END)
            self.entry_332.insert(0, "%.4f" % self.yclick)

        button332 = Button(self.frame2, text="Get",
                           command=getbutton332, highlightbackground='#2b2b2b', anchor=W, width=3)
        button332.grid(row=1, column=5, sticky=W)

        def CUT():
            self.lowercut = float(self.entry_31.get())
            self.highercut = float(self.entry_32.get())
            self.transcutlow = float(self.entry_33.get())
            self.transcut = float(self.entry_332.get())
            self.FTIRplot.set_xlim([self.lowercut, self.highercut])
            self.FTIRplot.set_ylim([self.transcutlow, self.transcut])
            self.canvas.show()

        def zoomall():
            self.entry_31.delete(0, END)
            self.entry_32.delete(0, END)
            self.entry_33.delete(0, END)
            self.entry_332.delete(0, END)
            self.entry_31.insert(0, 400)
            self.entry_32.insert(0, 6000)
            self.entry_33.insert(0, 0)
            self.entry_332.insert(0, 70)
            self.wavenumbers_cut = self.wavenumbers
            self.trans_cut = self.transmissions
            CUT()

        button34 = Button(self.frame2, text="Zoom all",
                          command=zoomall, highlightbackground='#2b2b2b', width=8)
        button34.grid(row=0, column=6)

        button35 = Button(self.frame2, text="CUT",
                          command=CUT, highlightbackground='#2b2b2b', width=8)
        button35.grid(row=1, column=6, sticky=W)

        label_21 = Label(self.frame2, text='Scale %', width=LABEL_WIDTH - 6, anchor=E, fg="#a9b7c6", bg='#2b2b2b')
        label_21.grid(row=0, column=7, columnspan=1, sticky=E)
        self.entry_21 = Entry(self.frame2, highlightbackground='#2b2b2b', width=4)
        self.entry_21.grid(row=0, column=8)
        self.entry_21.insert(0, "0.7")

        label_22 = Label(self.frame2, text='Angle:', width=LABEL_WIDTH - 6, anchor=E, fg="#a9b7c6", bg='#2b2b2b')
        label_22.grid(row=1, column=7, columnspan=1, sticky=E)
        self.entry_22 = Entry(self.frame2, highlightbackground='#2b2b2b', width=4)
        self.entry_22.grid(row=1, column=8)
        self.entry_22.insert(0, "0.0")

        label_23 = Label(self.frame2, text='CdTe offset(%):', width=LABEL_WIDTH, anchor=E, fg="#a9b7c6", bg='#2b2b2b')
        label_23.grid(row=0, column=9, columnspan=1, sticky=E)
        self.entry_23 = Entry(self.frame2, highlightbackground='#2b2b2b', width=4)
        self.entry_23.grid(row=0, column=10)
        self.entry_23.insert(0, "0.0")

        label_24 = Label(self.frame2, text='HgTe offset(%):', width=LABEL_WIDTH, anchor=E, fg="#a9b7c6", bg='#2b2b2b')
        label_24.grid(row=1, column=9, columnspan=1, sticky=E)
        self.entry_24 = Entry(self.frame2, highlightbackground='#2b2b2b', width=4)
        self.entry_24.grid(row=1, column=10)
        self.entry_24.insert(0, "0.0")

        self.intial_thicknesses_or_not = 1
        self.entry_d_list_initial = []

        button_21 = Button(self.frame2, text="Set",
                           command=self.setoffsets, highlightbackground='#2b2b2b', anchor=W, width=3)
        button_21.grid(row=1, column=11, sticky=W)

        self.frame1 = Frame(self, width=850, bg='#2b2b2b')
        self.frame1.pack(side=TOP, fill=BOTH, expand=True)
        # self.frame1.pack_propagate(0)

        self.FTIRfigure = Figure(figsize=(7.8, 4), dpi=100)
        self.FTIRfigure.subplots_adjust(left=0.08, bottom=0.12, right=0.92, top=0.95)
        self.FTIRplot = self.FTIRfigure.add_subplot(111)

        self.FTIRplot.plot(self.wavenumbers, self.transmissions)
        self.FTIRplot.set_xlim([self.lowercut, self.highercut])
        self.FTIRplot.set_ylim([0, self.transcut])
        self.FTIRplot.set_xlabel('Wavenumbers (cm-1)')
        self.FTIRplot.set_ylabel('Transmission (%)')
        self.FTIRplot.grid(True)
        self.vline = self.FTIRplot.axvline(x=400, visible=True, color='k', linewidth=0.7)
        self.hline = self.FTIRplot.axhline(y=0, visible=True, color='k', linewidth=0.7)
        self.dot = self.FTIRplot.plot(0, 0, marker='o', color='r')

        # a tk.DrawingArea
        self.canvas = FigureCanvasTkAgg(self.FTIRfigure, self.frame1)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.frame1)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        def on_key_event(event):
            # print('you pressed %s' % event.key)
            key_press_handler(event, self.canvas, self.toolbar)

        self.canvas.mpl_connect('key_press_event', on_key_event)
        self.canvas.mpl_connect('button_press_event', self.onpick)

        def load_structure_event(event):
            self.load_structure()

        def fit_fringes_event(event):
            self.fit_fringes()

        def openfromfile_event(event):
            self.openfromfile()

        def show_fringes_event(event):
            self.show_fringes()

        def cal_absorption_event(event):
            self.cal_absorption()

        def clearalldata_event(event):
            self.clearalldata()

        def help_event(event):
            self.help()

        if _platform == "darwin" or _platform == "linux" or _platform == "linux2":
            masterroot.bind('<Command-l>', load_structure_event)  # key must be binded to the tk window(unknown reason)
            masterroot.bind('<Command-g>', fit_fringes_event)
            masterroot.bind('<Command-o>', openfromfile_event)
            masterroot.bind('<Command-s>', show_fringes_event)
            masterroot.bind('<Command-a>', cal_absorption_event)
            masterroot.bind('<Command-c>', clearalldata_event)
            masterroot.bind('<Command-p>', help_event)
        elif _platform == "win32" or _platform == "win64":
            masterroot.bind('<Control-l>', load_structure_event)  # key must be binded to the tk window(unknown reason)
            masterroot.bind('<Control-g>', fit_fringes_event)
            masterroot.bind('<Control-o>', openfromfile_event)
            masterroot.bind('<Control-s>', show_fringes_event)
            masterroot.bind('<Control-a>', cal_absorption_event)
            masterroot.bind('<Control-c>', clearalldata_event)
            masterroot.bind('<Control-p>', help_event)

        self.pack()

    def help(self):
        self.addlog("FTIR fitting tool v. {}. With customization of layer structures. ".format(__version__))
        self.addlog('Open a FTIR transmission .csv file, then customize your layer structure on the right. '
                    'You can load or save a structure from file.')
        self.addlog('The "cut" and "zoom" function on the bottom can set the range of interest.')
        self.addlog('"scale factor" and "angle" settings are also on the bottom. Click "Show fringes" to see the result. ')
        self.addlog("The CdTe and HgTe offset function is used to fit the fringes. This is a very import idea in version 2. ")
        self.addlog("It is based on the fact that cell flux remains relatively stable during the growth. "
                    "So the layer thicknesses are not independent. They are related by the flux of CdTe and Te cells. ")
        self.addlog('Click "Set" to apply the offsets to the layer thicknesses. Note! Only the layers with check marks will be changed accordingly. ')
        self.addlog('Use"Fit fringes" to find the best CdTe/HgTe offset. Currently this function is running slow.')
        self.addlog('The fit range for CdTe offset is +-10, and for HgTe is +-5. You can change them inside the code. '
                    'The accuracy of the fitting is questionable due to fundamental reasons. ')
        self.addlog("A lot of functions currently are not working in version 2, check back on version 1 to see if they works. ")
        self.listbox.insert(END, '*' * 60)

    def setoffsets(self):

        """Adjust layer thicknesses using CdTe/HgTe offset. Only the layers with check marks will be changed. """

        if self.intial_thicknesses_or_not == 1:
            self.entry_d_list_initial = []

            for i in range(1, self.layernumber + 1):
                self.entry_d_list_initial.append(float(getattr(self, "entry_d_{}".format(i)).get()))
            self.intial_thicknesses_or_not = 0

        for i in range(1, self.layernumber + 1):
            if int(getattr(self, "checklayer{}".format(i)).get()) == 1:
                if getattr(self, "layertypevar{}".format(i)).get() == "CdTe":
                    new_d = self.entry_d_list_initial[i-1] * (1 + 0.01 * float(self.entry_23.get()))
                    getattr(self, "entry_d_{}".format(i)).delete(0, END)
                    getattr(self, "entry_d_{}".format(i)).insert(0, "{0:.2f}".format(new_d))
                elif getattr(self, "layertypevar{}".format(i)).get() == "MCT" or getattr(self, "layertypevar{}".format(i)).get() == "SL":
                    new_d = self.entry_d_list_initial[i-1] * float(getattr(self, "entry_x_{}".format(i)).get()) * (1 + 0.01 * float(self.entry_23.get())) \
                            + self.entry_d_list_initial[i-1] * (1 - float(getattr(self, "entry_x_{}".format(i)).get())) * (1 + 0.01 * float(self.entry_24.get()))
                    getattr(self, "entry_d_{}".format(i)).delete(0, END)
                    getattr(self, "entry_d_{}".format(i)).insert(0, "{0:.2f}".format(new_d))

    def settings(self):

        """Optinal settings for customized result."""

        settingwindow = Toplevel()
        w2 = 250  # width for the window
        h2 = 100  # height for the window
        ws = self.masterroot.winfo_screenwidth()  # width of the screen
        hs = self.masterroot.winfo_screenheight()  # height of the screen
        # calculate x and y coordinates for the Tk root window
        x2 = (ws / 2) - (w2 / 2)
        y2 = (hs / 3) - (h2 / 3)
        # set the dimensions of the screen
        # and where it is placed
        settingwindow.geometry('%dx%d+%d+%d' % (w2, h2, x2, y2))
        settingwindow.wm_title("Settings")
        # openfromfilewindow.wm_overrideredirect(True)
        settingwindow.configure(background='#2b2b2b', takefocus=True)
        settingwindow.attributes('-topmost', 'true')
        settingwindow.grab_set()

        Label(settingwindow, text="-------------Show Fringes------------", bg='#2b2b2b', fg="#a9b7c6", anchor=W).grid(row=0, column=0,
                                                                                             columnspan=2, sticky=W)

        self.displayreflection_temp, self.displayabsorption_temp = IntVar(), IntVar()
        checkboxr = Checkbutton(settingwindow, text="Show Reflection", variable=self.displayreflection_temp, bg='#2b2b2b', fg="#a9b7c6")
        checkboxr.grid(row=1, column=0, columnspan=2, sticky=W)

        if self.displayreflection == 1:
            checkboxr.select()

        checkboxa = Checkbutton(settingwindow, text="Show Absorption", variable=self.displayabsorption_temp, bg='#2b2b2b', fg="#a9b7c6")
        checkboxa.grid(row=2, column=0, columnspan=2, sticky=W)

        if self.displayabsorption== 1:
            checkboxa.select()

        # Samplenamegetoption.grab_set()
        # Structurenamegetoption.focus_set()

        def buttonOkayfuncton():
            self.displayreflection = self.displayreflection_temp.get()
            self.displayabsorption = self.displayabsorption_temp.get()

            settingwindow.grab_release()
            self.masterroot.focus_set()
            settingwindow.destroy()
            return

        def buttonCencelfuncton():
            settingwindow.grab_release()
            self.masterroot.focus_set()
            settingwindow.destroy()
            return

        def buttonOkayfunction_event(event):
            buttonOkayfuncton()

        buttonOK = Button(settingwindow, text="OK",
                          command=buttonOkayfuncton, highlightbackground='#2b2b2b', width=10)
        buttonOK.grid(row=3, column=0, columnspan=1)
        buttonOK = Button(settingwindow, text="Cancel",
                          command=buttonCencelfuncton, highlightbackground='#2b2b2b', width=10)
        buttonOK.grid(row=3, column=1, columnspan=1)
        settingwindow.bind('<Return>', buttonOkayfunction_event)

    def openfromfile(self):

        """Open a FTIR transmission .csv file. """

        if self.numberofdata >= 6:
            self.addlog('Cannot add more data file.')
            return

        self.wavenumbers = []
        self.transmissions = []
        self.absorptions = []

        file = filedialog.askopenfile(mode='r', defaultextension=".csv")
        if file is None:  # asksaveasfile return `None` if dialog closed with "cancel".
            return
        self.filename = file.name

        i = -1
        while self.filename[i] != '/':
            i -= 1
        self.filename = self.filename[i + 1:None]

        if self.filename[-4:None] != ".csv" and self.filename[-4:None] != ".CSV":
            self.addlog('{} format is not supported. Please select a .CSV file to open.'.format(self.filename[-4:None]))
            return

        self.filepath.config(text=self.filename)
        if self.filename[0:3] == 'Abs':
            with open(file.name, 'r') as f:
                reader = csv.reader(f, delimiter=',')
                for row in reader:
                    try:
                        self.wavenumbers.append(float(row[0]))
                        self.absorptions.append(float(row[1]))
                    except ValueError:
                        pass
            file.close()

            self.absorptionplot = self.FTIRplot.twinx()
            self.fitline_absorption = self.absorptionplot.plot(self.wavenumbers, self.absorptions, self.colororders2[self.numberofdata2], label=self.filename)
            self.absorptionplot.set_ylabel('Absorption Coefficient (cm-1)')
            self.absorptionplot.set_xlim([self.lowercut, self.highercut])
            self.absorptionplot.set_ylim([0, 10000])

            legend = self.absorptionplot.legend(loc='upper right', shadow=True)
            frame = legend.get_frame()
            frame.set_facecolor('0.90')

            # Set the fontsize
            for label in legend.get_texts():
                label.set_fontsize('medium')

            for label in legend.get_lines():
                label.set_linewidth(1.5)

            self.addlog('Added data {} ({})'.format(self.filename, self.colororders2[self.numberofdata2]))
            self.numberofdata2 += 1

        else:
            with open(file.name, 'r') as f:
                reader = csv.reader(f, delimiter=',')
                for row in reader:
                    try:
                        self.wavenumbers.append(float(row[0]))
                        self.transmissions.append(float(row[1]))
                    except ValueError:
                        pass
            file.close()

            # self.FTIRplot = self.FTIRfigure.add_subplot(111)
            self.FTIRplot.plot(self.wavenumbers, self.transmissions, self.colororders[self.numberofdata], label=self.filename)
            self.FTIRplot.set_xlim([self.lowercut, self.highercut])
            self.FTIRplot.set_ylim([self.transcutlow, self.transcut])
            self.FTIRplot.set_xlabel('Wavenumbers (cm-1)')
            self.FTIRplot.set_ylabel('Transmission (%)')
            self.FTIRplot.grid(True)

            legend = self.FTIRplot.legend(loc='upper right', shadow=True)
            frame = legend.get_frame()
            frame.set_facecolor('0.90')

            # Set the fontsize
            for label in legend.get_texts():
                label.set_fontsize('medium')

            for label in legend.get_lines():
                label.set_linewidth(1.5)

            # plt.savefig("test.png", dpi=300)
            # self.FTIRfigure.show()

            self.addlog('Added data {} ({})'.format(self.filename, self.colororders[self.numberofdata]))
            self.numberofdata += 1

        self.canvas.show()

        if len(self.wavenumbers) == 5810:
            self.addlog('Sample is probably characterized at EPIR.')
        elif len(self.wavenumbers) == 1946:
            self.addlog('Sample is probably characterized at UIC.')

    def savetofile(self):

        """Save calculated Transmission/Reflection/Absorption to file."""

        if self.peakvalues_fit ==[]:
            self.addlog("No data can be saved. ")
            return

        saveascsv = filedialog.asksaveasfilename(defaultextension='.csv')
        if saveascsv is None:
            return
        if saveascsv[-4:None] != ".csv" and saveascsv[-4:None] != ".CSV":
            self.addlog('Only .csv file can be saved.')
            return
        f = open(saveascsv, "w")
        if self.displayreflection == 1 or self.displayabsorption == 1:
            f.write("wn,T,R,A\n")

            for i in range(0, len(self.wavenumbers_cut)):
                f.write("{0:.6e},{1:.6e},{2:.6e},{3:.6e}\n".format(self.wavenumbers_cut[i], self.peakvalues_fit[i],
                                                                   self.reflections_fit[i], self.absorptions_fit[i]))
        else:
            f.write("wn,T\n")

            for i in range(0, len(self.wavenumbers_cut)):
                f.write("{0:.6e},{1:.6e}\n".format(self.wavenumbers_cut[i], self.peakvalues_fit[i]))

        f.close()

        self.addlog('Saved the file to: {}'.format(saveascsv))

    def load_structure(self):

        """Load existing heterojunction structures (.csv files.)"""

        self.osdir = "/Users/apple/Dropbox/6.python/Projects/FTIR_fittingtool_v2/Preload_Structure"

        filelist = []
        try:
            for item in os.listdir(self.osdir):
                if item[-4:None] == ".CSV" or item[-4:None] == ".csv":
                    filelist.append(item[0:-4])
        except FileNotFoundError:
            findornot = messagebox.askquestion(" ", "Structure folder not found. Do you want to relocate the folder manually?", icon='warning')
            if findornot == 'yes':
                self.osdir = filedialog.askdirectory()
                for item in os.listdir(self.osdir):
                    if item[-4:None] == ".CSV":
                        filelist.append(item[0:-4])
                self.addlog('Folder directory changed to {}'.format(self.osdir))
            else:
                return

        if filelist == []:
            self.addlog("No structure is found. ")
            self.osdir = "/Users/apple/Dropbox/6.python/Projects/FTIR_fittingtool_v2/Preload_Structure2"
            return

        openfromfilewindow = Toplevel()
        w2 = 250  # width for the window
        h2 = 80  # height for the window
        ws = self.masterroot.winfo_screenwidth()  # width of the screen
        hs = self.masterroot.winfo_screenheight()  # height of the screen
        # calculate x and y coordinates for the Tk root window
        x2 = (ws / 2) - (w2 / 2)
        y2 = (hs / 3) - (h2 / 3)
        # set the dimensions of the screen
        # and where it is placed
        openfromfilewindow.geometry('%dx%d+%d+%d' % (w2, h2, x2, y2))
        openfromfilewindow.wm_title("Open from file")
        # openfromfilewindow.wm_overrideredirect(True)
        openfromfilewindow.configure(background='#2b2b2b', takefocus=True)
        openfromfilewindow.attributes('-topmost', 'true')
        openfromfilewindow.grab_set()

        Label(openfromfilewindow, text="Existing Structure:", bg='#2b2b2b', fg="#a9b7c6", anchor=W).grid(row=0, column=0, columnspan=1, sticky=W)

        Structurenameget = StringVar(openfromfilewindow)
        Structurenameget.set("nBn_with_SL_barrier")  # initial value

        filelist = sorted(filelist)

        Structurenamegetoption = OptionMenu(openfromfilewindow, Structurenameget, *filelist)
        # '*'  to receive each list item as a separate parameter.
        Structurenamegetoption.config(bg='#2b2b2b')
        Structurenamegetoption.config(width=24)
        Structurenamegetoption["menu"].config(bg='#2b2b2b')
        Structurenamegetoption.grid(row=1, column=0, columnspan=2)
        # Samplenamegetoption.grab_set()
        Structurenamegetoption.focus_set()

        def buttonOpenfuncton():
            self.layernumber = 0

            self.Structurename = Structurenameget.get() + ".CSV"
            filename = self.osdir + "/" + self.Structurename

            layerlist = []
            xlist = []
            dlist = []
            checklist = []

            with open(filename, 'r') as f:
                reader = csv.reader(f, delimiter=',')
                for row in reader:
                    try:
                        layerlist.append(row[0])
                        xlist.append(float(row[1]))
                        dlist.append(float(row[2]))
                        checklist.append(int(row[3]))
                    except ValueError:
                        pass
            # file.close()

            def add_layer_on_top(layer, x, d, check_or_not):
                self.layernumber += 1

                self.buttonaddlayer.grid_forget()
                getattr(self, "layertypevar{}".format(self.layernumber)).set(layer)
                layertypeoption1 = OptionMenu(self.frame3, getattr(self, "layertypevar{}".format(self.layernumber)), *self.available_materials)
                layertypeoption1.config(bg="#2b2b2b", width=self.COLUMN0_WIDTH, anchor=E)
                layertypeoption1.grid(row=11 + self.layernumber, column=0, sticky=W+E)

                getattr(self, "entry_x_{}".format(self.layernumber)).grid(row=11 + self.layernumber, column=1, sticky=W+E)
                getattr(self, "entry_x_{}".format(self.layernumber)).insert(0, x)

                getattr(self, "entry_d_{}".format(self.layernumber)).grid(row=11 + self.layernumber, column=2)
                getattr(self, "entry_d_{}".format(self.layernumber)).insert(0, d)

                checkbox1 = Checkbutton(self.frame3, text="", variable=getattr(self, "checklayer{}".format(self.layernumber)), bg='#2b2b2b')
                checkbox1.grid(row=11 + self.layernumber, column=3, sticky=W)
                if check_or_not == 1:
                    checkbox1.select()

                if self.layernumber < 16:
                    self.buttonaddlayer.grid(row=12 + self.layernumber, column=0, columnspan=4, sticky=W+E)

            for i in range(0, len(layerlist)):
                if layerlist[i] in self.available_materials:
                    add_layer_on_top(layerlist[i], xlist[i], dlist[i], checklist[i])
                else:
                    self.addlog("Invalid Structure file.")
                    openfromfilewindow.grab_release()
                    self.masterroot.focus_set()
                    openfromfilewindow.destroy()
                    return

            self.addlog('Loaded Structure {}.'.format(self.Structurename[0:-4]))
            openfromfilewindow.grab_release()
            self.masterroot.focus_set()
            openfromfilewindow.destroy()
            return

        def buttonCencelfuncton():
            openfromfilewindow.grab_release()
            self.masterroot.focus_set()
            openfromfilewindow.destroy()
            return

        def buttonOpenfunction_event(event):
            buttonOpenfuncton()

        buttonOK = Button(openfromfilewindow, text="Open",
                           command=buttonOpenfuncton, highlightbackground='#2b2b2b', width=10)
        buttonOK.grid(row=2, column=0, columnspan=1)
        buttonOK = Button(openfromfilewindow, text="Cancel",
                           command=buttonCencelfuncton, highlightbackground='#2b2b2b', width=10)
        buttonOK.grid(row=2, column=1, columnspan=1)
        openfromfilewindow.bind('<Return>', buttonOpenfunction_event)

    def save_structure(self):

        """Save the customized structure to file. This is the best way to create structure files."""

        saveascsv = filedialog.asksaveasfilename(defaultextension='.CSV')
        if saveascsv is None:
            return
        if saveascsv[-4:None] != ".csv" and saveascsv[-4:None] != ".CSV":
            self.addlog('Only .csv file can be saved.')
            return
        f = open(saveascsv, "w")

        self.layertype_list, self.entry_x_list, self.entry_d_list, self.checklayer_list = [], [], [], []

        for i in range(1, self.layernumber + 1):
            self.layertype_list.append(getattr(self, "layertypevar{}".format(i)).get())
            self.entry_x_list.append(float(getattr(self, "entry_x_{}".format(i)).get()))
            self.entry_d_list.append(float(getattr(self, "entry_d_{}".format(i)).get()))
            self.checklayer_list.append(float(getattr(self, "checklayer{}".format(i)).get()))

        for i in range(0, len(self.layertype_list)):
            f.write("{},{},{},{}\n".format(self.layertype_list[i], self.entry_x_list[i], self.entry_d_list[i], int(self.checklayer_list[i])))

        f.close()

        self.addlog('Saved the Structure to: {}'.format(saveascsv))

    def show_fringes(self):

        """Show the calculated fringes curve based the structure and parameters provided. """

        self.layertype_list, self.entry_x_list, self.entry_d_list, self.checklayer_list = [], [], [], []

        for i in range(1, self.layernumber + 1):
            self.layertype_list.append(getattr(self, "layertypevar{}".format(i)).get())
            self.entry_x_list.append(float(getattr(self, "entry_x_{}".format(i)).get()))
            self.entry_d_list.append(float(getattr(self, "entry_d_{}".format(i)).get()))
            self.checklayer_list.append(int(getattr(self, "checklayer{}".format(i)).get()))

        self.wavenumbers_cut = []
        self.trans_cut = []

        addon = 0
        if self.wavenumbers == [] and self.transmissions == []:
            while addon <= 6000:
                self.wavenumbers_cut.append(500 + addon)
                self.trans_cut.append(0)
                addon += 5
        for i in range(0, len(self.wavenumbers)):
            if float(self.entry_32.get()) > float(self.wavenumbers[i]) > float(self.entry_31.get()):
                self.wavenumbers_cut.append(float(self.wavenumbers[i]))
                self.trans_cut.append(float(self.transmissions[i]))

        fitobject = FIT_FTIR(self.wavenumbers_cut, self.trans_cut, self.entry_d_0.get(), self.layertype_list, self.entry_x_list,
                             self.entry_d_list, self.checklayer_list, float(self.entry_21.get()), float(self.entry_22.get()),
                             float(self.entry_23.get()), float(self.entry_24.get()), self.subtype, 2,
                             self.listbox, self.progress_var, self.wn_beingcalculated)
        self.peakvalues_fit = fitobject.returnpeakvalues()
        self.reflections_fit = fitobject.returnreflections()
        self.absorptions_fit = fitobject.returnabsorptions()
        
        if self.displayreflection == 0 and self.displayabsorption == 0:
            try:
                self.fitline2.pop(0).remove()
                self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'r')
            except (AttributeError, IndexError) as error:
                try:
                    self.fitline.pop(0).remove()
                    self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'r')
                except (AttributeError, IndexError) as error:
                    self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'r')

            try:
                self.fitline3.pop(0).remove()
            except (AttributeError, IndexError) as error:
                pass

            try:
                self.fitline4.pop(0).remove()
            except (AttributeError, IndexError) as error:
                pass

            self.FTIRplot.set_xlim([self.lowercut, self.highercut])
            self.FTIRplot.set_ylim([self.transcutlow, self.transcut])
            self.FTIRplot.set_xlabel('Wavenumbers (cm-1)')
            self.FTIRplot.set_ylabel('Transmission (%)')
            self.FTIRplot.grid(True)
            self.canvas.show()
        elif self.displayreflection == 1 and self.displayabsorption == 0:
            try:
                self.fitline2.pop(0).remove()
                self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'g')
            except (AttributeError, IndexError) as error:
                try:
                    self.fitline.pop(0).remove()
                    self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'g')
                except (AttributeError, IndexError) as error:
                    self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'g')

            try:
                self.fitline3.pop(0).remove()
                self.fitline3 = self.FTIRplot.plot(self.wavenumbers_cut, self.reflections_fit, 'r')
            except (AttributeError, IndexError) as error:
                try:
                    self.fitline.pop(0).remove()
                    self.fitline3 = self.FTIRplot.plot(self.wavenumbers_cut, self.reflections_fit, 'r')
                except (AttributeError, IndexError) as error:
                    self.fitline3 = self.FTIRplot.plot(self.wavenumbers_cut, self.reflections_fit, 'r')

            try:
                self.fitline4.pop(0).remove()
            except (AttributeError, IndexError) as error:
                pass

            self.FTIRplot.set_xlim([self.lowercut, self.highercut])
            self.FTIRplot.set_ylim([self.transcutlow, self.transcut])
            self.FTIRplot.set_xlabel('Wavenumbers (cm-1)')
            self.FTIRplot.set_ylabel('Transmission/Reflection (%)')
            self.FTIRplot.grid(True)
            self.canvas.show()

        elif self.displayreflection == 0 and self.displayabsorption == 1:
            try:
                self.fitline2.pop(0).remove()
                self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'g')
            except (AttributeError, IndexError) as error:
                try:
                    self.fitline.pop(0).remove()
                    self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'g')
                except (AttributeError, IndexError) as error:
                    self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'g')

            try:
                self.fitline4.pop(0).remove()
                self.fitline4 = self.FTIRplot.plot(self.wavenumbers_cut, self.absorptions_fit, 'purple')
            except (AttributeError, IndexError) as error:
                try:
                    self.fitline.pop(0).remove()
                    self.fitline4 = self.FTIRplot.plot(self.wavenumbers_cut, self.absorptions_fit, 'purple')
                except (AttributeError, IndexError) as error:
                    self.fitline4 = self.FTIRplot.plot(self.wavenumbers_cut, self.absorptions_fit, 'purple')

            self.FTIRplot.set_xlim([self.lowercut, self.highercut])
            self.FTIRplot.set_ylim([self.transcutlow, self.transcut])
            self.FTIRplot.set_xlabel('Wavenumbers (cm-1)')
            self.FTIRplot.set_ylabel('Transmission/Absorption (%)')
            self.FTIRplot.grid(True)
            self.canvas.show()

            try:
                self.fitline3.pop(0).remove()
            except (AttributeError, IndexError) as error:
                pass

        elif self.displayreflection == 1 and self.displayabsorption == 1:
            try:
                self.fitline2.pop(0).remove()
                self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'g')
            except (AttributeError, IndexError) as error:
                try:
                    self.fitline.pop(0).remove()
                    self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'g')
                except (AttributeError, IndexError) as error:
                    self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'g')

            try:
                self.fitline3.pop(0).remove()
                self.fitline3 = self.FTIRplot.plot(self.wavenumbers_cut, self.reflections_fit, 'r')
            except (AttributeError, IndexError) as error:
                try:
                    self.fitline.pop(0).remove()
                    self.fitline3 = self.FTIRplot.plot(self.wavenumbers_cut, self.reflections_fit, 'r')
                except (AttributeError, IndexError) as error:
                    self.fitline3 = self.FTIRplot.plot(self.wavenumbers_cut, self.reflections_fit, 'r')

            try:
                self.fitline4.pop(0).remove()
                self.fitline4 = self.FTIRplot.plot(self.wavenumbers_cut, self.absorptions_fit, 'purple')
            except (AttributeError, IndexError) as error:
                try:
                    self.fitline.pop(0).remove()
                    self.fitline4 = self.FTIRplot.plot(self.wavenumbers_cut, self.absorptions_fit, 'purple')
                except (AttributeError, IndexError) as error:
                    self.fitline4 = self.FTIRplot.plot(self.wavenumbers_cut, self.absorptions_fit, 'purple')

            self.FTIRplot.set_xlim([self.lowercut, self.highercut])
            self.FTIRplot.set_ylim([self.transcutlow, self.transcut])
            self.FTIRplot.set_xlabel('Wavenumbers (cm-1)')
            self.FTIRplot.set_ylabel('Transmission/Reflection/Absorption (%)')
            self.FTIRplot.grid(True)
            self.canvas.show()

    def fit_fringes(self):

        """Find the best CdTe/HgTe offsets to fit the fringes. """

        self.layertype_list, self.entry_x_list, self.entry_d_list, self.checklayer_list = [], [], [], []

        for i in range(1, self.layernumber + 1):
            self.layertype_list.append(getattr(self, "layertypevar{}".format(i)).get())
            self.entry_x_list.append(float(getattr(self, "entry_x_{}".format(i)).get()))
            self.entry_d_list.append(float(getattr(self, "entry_d_{}".format(i)).get()))
            self.checklayer_list.append(int(getattr(self, "checklayer{}".format(i)).get()))

        self.wavenumbers_cut = []
        self.trans_cut = []

        if float(self.entry_32.get()) > 5000:
            self.addlog('Please choose the cut range of fringes.')
            return

        for i in range(0, len(self.wavenumbers)):
            if float(self.entry_32.get()) > float(self.wavenumbers[i]) > float(self.entry_31.get()):
                self.wavenumbers_cut.append(float(self.wavenumbers[i]))
                self.trans_cut.append(float(self.transmissions[i]))

        if self.intial_thicknesses_or_not == 1:
            self.entry_d_list_initial = []

            for i in range(1, self.layernumber + 1):
                self.entry_d_list_initial.append(float(getattr(self, "entry_d_{}".format(i)).get()))
            self.intial_thicknesses_or_not = 0

        self.addprogressbar()
        self.text2 = self.status2.cget("text")

        self.addlog('*' * 60)
        self.addlog("Fitting fringes in process. Please wait...")

        self.queue = queue.Queue()

        ThreadedTask_fringes(self.queue, float(self.entry_23.get()), float(self.entry_24.get()), self.entry_d_list_initial,
                      self.layernumber, self.wavenumbers_cut, self.trans_cut, float(self.entry_d_0.get()),
                      self.layertype_list, self.entry_x_list, self.entry_d_list, self.checklayer_list,
                      float(self.entry_21.get()), float(self.entry_22.get()),
                      self.subtype, 2, self.listbox, self.progress_var, self.wn_beingcalculated).start()
        self.master.after(100, self.process_queue_fringes)

    def process_queue_fringes(self):

        """Threading function for self.fit_fringes()."""

        try:
            result = self.queue.get(0)
            # Show result of the task if needed
            self.entry_23.delete(0, END)
            self.entry_23.insert(0, '{0:.2f}'.format(result[0]))
            self.entry_24.delete(0, END)
            self.entry_24.insert(0, '{0:.2f}'.format(result[1]))

            self.setoffsets()

            self.layertype_list, self.entry_x_list, self.entry_d_list, self.checklayer_list = [], [], [], []

            for i in range(1, self.layernumber + 1):
                self.layertype_list.append(getattr(self, "layertypevar{}".format(i)).get())
                self.entry_x_list.append(float(getattr(self, "entry_x_{}".format(i)).get()))
                self.entry_d_list.append(float(getattr(self, "entry_d_{}".format(i)).get()))
                self.checklayer_list.append(int(getattr(self, "checklayer{}".format(i)).get()))

            fitobject = FIT_FTIR(self.wavenumbers_cut, self.trans_cut, self.entry_d_0.get(), self.layertype_list, self.entry_x_list,
                                 self.entry_d_list, self.checklayer_list, float(self.entry_21.get()), float(self.entry_22.get()),
                                 float(self.entry_23.get()), float(self.entry_24.get()), self.subtype, 2,
                                 self.listbox, self.progress_var, self.wn_beingcalculated)
            self.peakvalues_fit = fitobject.returnpeakvalues()

            self.addlog('Fitting fringes complete. MSE={}'.format(result[2]))

            try:
                self.fitline2.pop(0).remove()
                self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'r')
            except (AttributeError, IndexError) as error:
                try:
                    self.fitline.pop(0).remove()
                    self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'r')
                except (AttributeError, IndexError) as error:
                    self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'r')

            self.FTIRplot.set_xlim([self.lowercut, self.highercut])
            self.FTIRplot.set_ylim([self.transcutlow, self.transcut])
            self.FTIRplot.set_xlabel('Wavenumbers (cm-1)')
            self.FTIRplot.set_ylabel('Transmission (%)')
            self.FTIRplot.grid(True)
            self.canvas.show()

            self.removeprogressbar()

        except queue.Empty:
            self.after(100, self.process_queue_fringes)

    def cal_absorption(self):

        """Calculate the absorption coefficient as a function of wavenumebrs. """

        self.wavenumbers_cut1 = []
        self.trans_cut1 = []
        skipnumber = int(len(self.wavenumbers) / 200)

        for i in range(0, len(self.wavenumbers)):
            if float(self.entry_32.get()) >= float(self.wavenumbers[i]) >= float(self.entry_31.get()):
                self.wavenumbers_cut1.append(float(self.wavenumbers[i]))
                if self.transmissions[i] >= 0:
                    self.trans_cut1.append(float(self.transmissions[i]))
                else:
                    self.trans_cut1.append(0)

        roughornot = messagebox.askquestion(" ", "Do a rough calculation?", icon='warning')
        if roughornot == 'yes':
            wavenumbers_cut_rough = []
            trans_cut_rough = []
            index1 = 0
            while index1 in range(0, len(self.wavenumbers_cut1)):
                wavenumbers_cut_rough.append(self.wavenumbers_cut1[index1])
                trans_cut_rough.append(self.trans_cut1[index1])
                index1 += skipnumber
            self.wavenumbers_cut1 = wavenumbers_cut_rough
            self.trans_cut1 = trans_cut_rough
        else:
            pass

        self.addprogressbar()
        self.text2 = self.status2.cget("text")

        self.addlog('*' * 60)
        self.addlog("Absorption calculation in process. Please wait...")

        self.queue = queue.Queue()
        ThreadedTask1(self.queue, float(self.entry_24.get()), float(self.entry_29.get()), self.wavenumbers_cut1,
                      self.trans_cut1, float(self.entry_25.get()), float(self.entry_26.get()),
                      float(self.entry_27.get()),
                      float(self.entry_28.get()), float(self.entry_290.get()), float(self.entry_291.get()),
                      float(self.entry_292.get()), float(self.entry_293.get()), float(self.entry_294.get()),
                      float(self.entry_295.get()), float(self.entry_296.get()), float(self.entry_21.get()),
                      self.subtype, 0, self.listbox, self.progress_var, self.wn_beingcalculated).start()
        self.master.after(100, self.process_queue_absorption)

    def process_queue_absorption(self):

        """Threading function for self.cal_absorption()."""

        try:
            self.trackwavenumber()
            result = self.queue.get(0)
            # Show result of the task if needed
            self.absorptions = result

            self.addlog('Absorption calculation complete!')

            try:
                self.fitline2.pop(0).remove()
            except (AttributeError, IndexError) as error:
                pass

            self.absorptionplot = self.FTIRplot.twinx()
            self.fitline_absorption = self.absorptionplot.plot(self.wavenumbers_cut1, self.absorptions, self.colororders2[self.numberofdata2], label='Calculated Absorption')
            self.absorptionplot.set_ylabel('Absorption Coefficient (cm-1)')
            self.absorptionplot.set_xlim([self.lowercut, self.highercut])
            self.absorptionplot.set_ylim([0, 10000])

            legend = self.absorptionplot.legend(loc='upper right', shadow=True)
            frame = legend.get_frame()
            frame.set_facecolor('0.90')

            # Set the fontsize
            for label in legend.get_texts():
                label.set_fontsize('medium')

            for label in legend.get_lines():
                label.set_linewidth(1.5)

            self.canvas.show()

            self.numberofdata2 += 1

            saveornot = messagebox.askquestion(" ", "Save the result as a .csv file?", icon='warning')
            if saveornot == 'yes':
                saveascsv = filedialog.asksaveasfilename(defaultextension='.csv')
                if saveascsv is None:
                    self.removeprogressbar()
                    self.removewavenumber()
                    return
                if saveascsv[-4:None] != ".csv" and saveascsv[-4:None] != ".CSV":
                    self.addlog('Only .csv file can be saved.')
                    self.removeprogressbar()
                    self.removewavenumber()
                    return
                f = open(saveascsv, "w")

                for i in range(0, len(self.wavenumbers_cut1)):
                    f.write("{0:.6e},{1:.6e}\n".format(self.wavenumbers_cut1[i], self.absorptions[i]))

                f.close()

                self.addlog('Saved the file to: {}'.format(saveascsv))
                self.removeprogressbar()
                self.removewavenumber()
            else:
                self.removeprogressbar()
                self.removewavenumber()
                return
        except queue.Empty:
            self.after(100, self.process_queue_absorption)

    def cal_MCT_absorption(self):
        call_MCT_choosewindow = Toplevel()
        w2 = 250  # width for the window
        h2 = 80  # height for the window
        ws = self.masterroot.winfo_screenwidth()  # width of the screen
        hs = self.masterroot.winfo_screenheight()  # height of the screen
        # calculate x and y coordinates for the Tk root window
        x2 = (ws / 2) - (w2 / 2)
        y2 = (hs / 3) - (h2 / 3)
        # set the dimensions of the screen
        # and where it is placed
        call_MCT_choosewindow.geometry('%dx%d+%d+%d' % (w2, h2, x2, y2))
        call_MCT_choosewindow.wm_title("Choose fitting method")
        # call_MCT_choosewindow.wm_overrideredirect(True)
        call_MCT_choosewindow.configure(background='#2b2b2b', takefocus=True)
        call_MCT_choosewindow.attributes('-topmost', 'true')
        call_MCT_choosewindow.grab_set()

        Label(call_MCT_choosewindow, text="Use fitting method by:", bg='#2b2b2b', fg="#a9b7c6", anchor=W).grid(row=0, column=0, columnspan=1, sticky=W)

        methodnameget = StringVar(call_MCT_choosewindow)
        methodnameget.set("Yong")  # initial value

        methodnamegetoption = OptionMenu(call_MCT_choosewindow, methodnameget, "Chu", "Schacham and Finkman", "Yong")
        # '*'  to receive each list item as a separate parameter.
        methodnamegetoption.config(bg='#2b2b2b')
        methodnamegetoption.config(width=24)
        methodnamegetoption["menu"].config(bg='#2b2b2b')
        methodnamegetoption.grid(row=1, column=0, columnspan=2)
        # Samplenamegetoption.grab_set()
        methodnamegetoption.focus_set()

        def buttongofuncton():
            fitobject = cal_MCT_a(self.entry_24.get(), self.wavenumbers, methodnameget.get())

            self.absorptions = fitobject.return_absorptions()

            try:
                self.fitline_MCT.pop(0).remove()
            except (AttributeError, IndexError) as error:
                pass

            self.absorptionplot = self.FTIRplot.twinx()
            self.fitline_MCT = self.absorptionplot.plot(self.wavenumbers, self.absorptions,
                                                            self.colororders2[self.numberofdata2])
            self.absorptionplot.set_ylabel('Absorption Coefficient (cm-1)')
            self.absorptionplot.set_xlim([self.lowercut, self.highercut])
            self.absorptionplot.set_ylim([0, 10000])

            self.canvas.show()

            self.numberofdata2 += 1

            self.addlog("Showing MCT abosorption curve "
                        "for x = {} using {}'s formula(1994).".format(self.entry_24.get(), methodnameget.get()))

            call_MCT_choosewindow.grab_release()
            self.masterroot.focus_set()
            call_MCT_choosewindow.destroy()

            saveornot = messagebox.askquestion(" ", "Save the result as a .csv file?", icon='warning')
            if saveornot == 'yes':
                saveascsv = filedialog.asksaveasfilename(defaultextension='.csv')
                if saveascsv is None:
                    return
                if saveascsv[-4:None] != ".csv" and saveascsv[-4:None] != ".CSV":
                    self.addlog('Only .csv file can be saved.')
                    return
                f = open(saveascsv, "w")

                for i in range(0, len(self.wavenumbers)):
                    f.write("{0:.6e},{1:.6e}\n".format(self.wavenumbers[i], self.absorptions[i]))

                f.close()

                self.addlog('Saved the file to: {}'.format(saveascsv))
            else:
                return

        def buttonCancelfuncton():
            call_MCT_choosewindow.grab_release()
            self.masterroot.focus_set()
            call_MCT_choosewindow.destroy()
            return

        def buttongofunction_event(event):
            buttongofuncton()

        buttonOK = Button(call_MCT_choosewindow, text="Calculate",
                          command=buttongofuncton, highlightbackground='#2b2b2b', width=10)
        buttonOK.grid(row=2, column=0, columnspan=1, sticky=W)

        buttonOK = Button(call_MCT_choosewindow, text="Cancel",
                          command=buttonCancelfuncton, highlightbackground='#2b2b2b', width=10)
        buttonOK.grid(row=2, column=1, columnspan=1, sticky=E)
        call_MCT_choosewindow.bind('<Return>', buttongofunction_event)

    def fit_fringes_old(self):
        self.layertype_list, self.entry_x_list, self.entry_d_list, self.checklayer_list = [], [], [], []

        for i in range(1, self.layernumber + 1):
            self.layertype_list.append(getattr(self, "layertypevar{}".format(i)).get())
            self.entry_x_list.append(float(getattr(self, "entry_x_{}".format(i)).get()))
            self.entry_d_list.append(float(getattr(self, "entry_d_{}".format(i)).get()))
            self.checklayer_list.append(int(getattr(self, "checklayer{}".format(i)).get()))

        self.wavenumbers_cut = []
        self.trans_cut = []

        if float(self.entry_32.get()) > 5000:
            self.addlog('Please choose the cut range of fringes.')
            return

        for i in range(0, len(self.wavenumbers)):
            if float(self.entry_32.get()) > float(self.wavenumbers[i]) > float(self.entry_31.get()):
                self.wavenumbers_cut.append(float(self.wavenumbers[i]))
                self.trans_cut.append(float(self.transmissions[i]))

        def fit_func(x, CdTe_offset, HgTe_offset):
            self.entry_23.delete(0, END)
            self.entry_23.insert(0, '{0:.2f}'.format(CdTe_offset))
            self.entry_24.delete(0, END)
            self.entry_24.insert(0, '{0:.2f}'.format(HgTe_offset))

            self.setoffsets()

            self.layertype_list, self.entry_x_list, self.entry_d_list, self.checklayer_list = [], [], [], []

            for i in range(1, self.layernumber + 1):
                self.layertype_list.append(getattr(self, "layertypevar{}".format(i)).get())
                self.entry_x_list.append(float(getattr(self, "entry_x_{}".format(i)).get()))
                self.entry_d_list.append(float(getattr(self, "entry_d_{}".format(i)).get()))
                self.checklayer_list.append(int(getattr(self, "checklayer{}".format(i)).get()))

            lamda = 10000 / x
            fitobject1 = FIT_FTIR(self.wavenumbers_cut, self.trans_cut, self.entry_d_0.get(), self.layertype_list, self.entry_x_list,
                             self.entry_d_list, self.checklayer_list, float(self.entry_21.get()), float(self.entry_22.get()),
                             CdTe_offset, HgTe_offset, self.subtype, 0,
                             self.listbox, self.progress_var, self.wn_beingcalculated)
            return fitobject1.cal_fringes(lamda)

        p0 = np.array([float(self.entry_23.get()), float(self.entry_24.get())])

        self.popt, self.pcov = curve_fit(fit_func, self.wavenumbers_cut, self.trans_cut, p0,
                                         bounds=([float(self.entry_23.get()) - 5, float(self.entry_24.get()) - 5],
                                                 [float(self.entry_23.get()) + 5, float(self.entry_24.get()) + 5]))

        self.entry_23.delete(0, END)
        self.entry_23.insert(0, '{0:.2f}'.format(self.popt[0]))
        self.entry_24.delete(0, END)
        self.entry_24.insert(0, '{0:.2f}'.format(self.popt[1]))

        self.setoffsets()

        self.layertype_list, self.entry_x_list, self.entry_d_list, self.checklayer_list = [], [], [], []

        for i in range(1, self.layernumber + 1):
            self.layertype_list.append(getattr(self, "layertypevar{}".format(i)).get())
            self.entry_x_list.append(float(getattr(self, "entry_x_{}".format(i)).get()))
            self.entry_d_list.append(float(getattr(self, "entry_d_{}".format(i)).get()))
            self.checklayer_list.append(int(getattr(self, "checklayer{}".format(i)).get()))

        fitobject = FIT_FTIR(self.wavenumbers_cut, self.trans_cut, self.entry_d_0.get(), self.layertype_list, self.entry_x_list,
                             self.entry_d_list, self.checklayer_list, float(self.entry_21.get()), float(self.entry_22.get()),
                             float(self.entry_23.get()), float(self.entry_24.get()), self.subtype, 2,
                             self.listbox, self.progress_var, self.wn_beingcalculated)
        self.peakvalues_fit = fitobject.returnpeakvalues()

        self.addlog('Fitting fringes complete.')

        try:
            self.fitline2.pop(0).remove()
            self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'r')
        except (AttributeError, IndexError) as error:
            try:
                self.fitline.pop(0).remove()
                self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'r')
            except (AttributeError, IndexError) as error:
                self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'r')

        self.FTIRplot.set_xlim([self.lowercut, self.highercut])
        self.FTIRplot.set_ylim([0, self.transcut])
        self.FTIRplot.set_xlabel('Wavenumbers (cm-1)')
        self.FTIRplot.set_ylabel('Transmission (%)')
        self.FTIRplot.grid(True)
        self.canvas.show()

    def fit_fringes_old(self):
        self.wavenumbers_cut = []
        self.MSE = 10000000000
        self.MSE_new = 0
        self.fittedthickness = 0
        if float(self.entry_32.get()) > 5000:
            self.addlog('Please choose the cut range of fringes.')
            return

        for wn in self.wavenumbers:
            if float(self.entry_32.get()) > float(wn) > float(self.entry_31.get()):
                self.wavenumbers_cut.append(wn)
        for layerthickness in np.arange(float(self.entry_29.get()) - 4, float(self.entry_29.get()) + 4, 0.01):
            fitobject = FIT_FTIR(float(self.entry_24.get()), layerthickness, self.wavenumbers_cut,
                                 float(self.entry_25.get()), float(self.entry_26.get()), float(self.entry_27.get()),
                                 float(self.entry_28.get()), float(self.entry_290.get()), float(self.entry_291.get()),
                                 float(self.entry_292.get()), float(self.entry_293.get()), float(self.entry_294.get()),
                                 float(self.entry_21.get()), self.subtype, 2, self.listbox, self.progress_var, self.wn_beingcalculated)
            self.peakvalues_fit = fitobject.returnpeakvalues()

            for i in range(0, len(self.wavenumbers_cut)):
                self.MSE_new += 1 / len(self.wavenumbers_cut) * (
                float(self.peakvalues_fit[i]) - float(self.transmissions[i])) \
                                * (float(self.peakvalues_fit[i]) - float(self.transmissions[i]))

            if self.MSE_new < self.MSE:
                self.MSE = self.MSE_new
                self.fittedthickness = layerthickness
            self.MSE_new = 0

        self.entry_29.delete(0, END)
        self.entry_29.insert(0, "%.2f" % self.fittedthickness)

        fitobject = FIT_FTIR(float(self.entry_24.get()),
                             self.fittedthickness, self.wavenumbers_cut, float(self.entry_28.get()),
                             float(self.entry_290.get()), float(self.entry_291.get()), float(self.entry_292.get()),
                             float(self.entry_293.get()), float(self.entry_21.get()), self.subtype, 2, self.listbox,
                             self.progress_var, self.wn_beingcalculated)
        self.peakvalues_fit = fitobject.returnpeakvalues()

        self.addlog('Fitting fringes complete. MSE={}'.format(self.MSE))

        try:
            self.fitline2.pop(0).remove()
            self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'r')
        except (AttributeError, IndexError) as error:
            try:
                self.fitline.pop(0).remove()
                self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'r')
            except (AttributeError, IndexError) as error:
                self.fitline2 = self.FTIRplot.plot(self.wavenumbers_cut, self.peakvalues_fit, 'r')

        self.FTIRplot.set_xlim([self.lowercut, self.highercut])
        self.FTIRplot.set_ylim([0, self.transcut])
        self.FTIRplot.set_xlabel('Wavenumbers (cm-1)')
        self.FTIRplot.set_ylabel('Transmission (%)')
        self.FTIRplot.grid(True)
        self.canvas.show()

    def onpick(self, event):

        """Define the behaviors when the mouse click on the graph. """

        self.xclick = event.xdata
        self.yclick = event.ydata
        self.vline.set_xdata(self.xclick)
        self.hline.set_ydata(self.yclick)
        try:
            self.dot.pop(0).remove()
        except IndexError:
            if self.xclick is not None and self.yclick is not None:
                self.dot = self.FTIRplot.plot(self.xclick, self.yclick, marker='x', color='r')
            return
        if self.xclick is not None and self.yclick is not None:
            self.dot = self.FTIRplot.plot(self.xclick, self.yclick, marker='x', color='r')

        self.canvas.draw()

        def find_x():  # Find x knowing cutoff energy E
            self.composition = 0
            Delta = 100
            while Delta > 0.0001:
                Delta = self.energy - (
                    (-.302 + 1.93 * self.composition +
                     5.35 * .0001 * self.Temp * (1 - 2 * self.composition) -
                     .81 * self.composition * self.composition +
                     .832 * self.composition * self.composition * self.composition) * 1000)
                self.composition += 0.00001

        if self.xclick is not None and self.yclick is not None:
            self.wavenumber = self.xclick
            self.trans = self.yclick
            self.wavelength = 10000 / self.wavenumber
            self.energy = 4.13566743 * 3 * 100 / self.wavelength
            find_x()
        else:
            self.trans = 0
            self.wavenumber = 0
            self.wavelength = 0
            self.energy = 0
            self.composition = 0

        self.wavenumber1.config(text='{0:.4f}'.format(self.wavenumber))
        self.wavelength1.config(text='{0:.4f}'.format(self.wavelength))
        self.energy1.config(text='{0:.4f}'.format(self.energy))
        self.composition1.config(text='{0:.4f}'.format(self.composition))

        fitobject = FIT_FTIR(float(self.entry_24.get()), float(self.entry_29.get()), self.wavenumbers,
                             self.transmissions,
                             float(self.entry_25.get()), float(self.entry_26.get()), float(self.entry_27.get()),
                             float(self.entry_28.get()), float(self.entry_290.get()), float(self.entry_291.get()),
                             float(self.entry_292.get()), float(self.entry_293.get()), float(self.entry_294.get()),
                             float(self.entry_295.get()), float(self.entry_296.get()), float(self.entry_21.get()),
                             self.subtype, 0, self.listbox, self.progress_var, self.wn_beingcalculated)
        try:
            self.addlog('Absorption Coefficient: {}cm-1'.format(fitobject.cal_absorption_single(self.xclick)))
        except IndexError:
            pass

    def clearalldata(self):

        """Clear everything. """

        if self.wavenumber != 0 or self.filepath.cget('text') != '':
            self.FTIRplot.clear()
            self.absorptionplot.clear()
            self.lowercut = 400
            self.highercut = 6000
            self.transcut = 70
            self.FTIRplot.set_xlim([self.lowercut, self.highercut])
            self.FTIRplot.set_ylim([0, self.transcut])
            self.FTIRplot.set_xlabel('Wavenumbers (cm-1)')
            self.FTIRplot.set_ylabel('Transmission (%)')
            self.FTIRplot.grid(True)
            self.dot = self.FTIRplot.plot(0, 0, marker='x', color='w')
            self.vline = self.FTIRplot.axvline(x=400, visible=True, color='k', linewidth=0.7)
            self.hline = self.FTIRplot.axhline(y=0, visible=True, color='k', linewidth=0.7)
            self.canvas.show()
            self.filepath.config(text='')
            self.filename = ''

            self.trans = 0
            self.wavenumber = 0
            self.wavelength = 0
            self.energy = 0
            self.composition = 0
            self.wavenumber1.config(text='{}'.format(self.wavenumber))
            self.wavelength1.config(text='{}'.format(self.wavelength))
            self.energy1.config(text='{}'.format(self.energy))
            self.composition1.config(text='{}'.format(self.composition))

            self.entry_21.delete(0, END)
            self.entry_22.delete(0, END)
            self.entry_23.delete(0, END)
            self.entry_24.delete(0, END)
            self.entry_25.delete(0, END)
            self.entry_26.delete(0, END)
            self.entry_27.delete(0, END)
            self.entry_28.delete(0, END)
            self.entry_29.delete(0, END)
            self.entry_290.delete(0, END)
            self.entry_291.delete(0, END)
            self.entry_31.delete(0, END)
            self.entry_32.delete(0, END)
            self.entry_33.delete(0, END)
            self.entry_21.insert(0, '0.00')
            self.entry_22.insert(0, '500')
            self.entry_23.insert(0, '0.00')
            self.entry_24.insert(0, '0.00')
            self.entry_25.insert(0, '0.00')
            self.entry_26.insert(0, '0.00')
            self.entry_27.insert(0, '0.00')
            self.entry_28.insert(0, '0.00')
            self.entry_29.insert(0, '0.00')
            self.entry_290.insert(0, '10.00')
            self.entry_291.insert(0, '1.00')
            self.entry_31.insert(0, 400)
            self.entry_32.insert(0, 6000)
            self.entry_33.insert(0, 70)

            self.numberofdata = 0
            self.numberofdata2 = 0
            self.addlog('*' * 60)

    def addlog(self, string):
        self.listbox.insert(END, string)
        self.listbox.yview(END)

    def addprogressbar(self):
        self.text = self.status1.cget("text")
        self.status1.pack_forget()

        self.progressbar = Progressbar(self.statusbar, variable=self.progress_var, maximum=100)
        self.progressbar.pack(side=LEFT, fill=X, expand=1)

    def removeprogressbar(self):
        self.progressbar.pack_forget()

        self.status1 = Label(self.statusbar, text=self.text, fg='#a9b7c6', bg='#2b2b2b', bd=1, relief=RIDGE)
        self.status1.pack(side=LEFT, fill=X, expand=True)
        self.status1.pack_propagate(0)

    def trackwavenumber(self):
        self.status2.config(text='Wavenumber = {:.1f}cm-1'.format(self.wn_beingcalculated.get()))

    def removewavenumber(self):
        self.status2.config(text=self.text2)


def main():
    root = Tk()
    w = cross_platform_config.config.FRAME_WIDTH  # width for the Tk root
    h = cross_platform_config.config.FRAME_HEIGHT  # height for the Tk root
    ws = root.winfo_screenwidth()  # width of the screen
    hs = root.winfo_screenheight()  # height of the screen
    x = (ws / 2) - (w / 2)
    y = (hs / 4) - (h / 4)
    root.geometry('%dx%d+%d+%d' % (w, h, x, y))

    root.wm_title("FTIR Fitting Tool v. {}".format(__version__))
    # root.iconbitmap("icon.icns")
    root.configure(background='#2b2b2b')

    # Status bar #
    statusbar = Frame(root, bg='#2b2b2b', bd=1, relief=RIDGE)

    authorLabel = Label(statusbar, text='© 1,2018 Peihong Man.', fg='#a9b7c6', bg='#2b2b2b', bd=1, relief=RIDGE,
                        padx=4.2, width=21)
    authorLabel.pack(side=LEFT)
    authorLabel.pack_propagate(0)

    status1 = Label(statusbar, text='Welcome to FTIR Fitting Tool. Press ⌘+P for help.', fg='#a9b7c6', bg='#2b2b2b', bd=1, relief=RIDGE)
    if _platform == "win32" or _platform == "win64":
        status1.config(text='Welcome to FTIR Fitting Tool. Press Ctrl + P for help.')
    status1.pack(side=LEFT, fill=X, expand=True)
    status1.pack_propagate(0)
    status2 = Label(statusbar, text='v. {}'.format(__version__), fg='#a9b7c6', bg='#2b2b2b', bd=1, relief=RIDGE,
                    width=21)
    status2.pack(side=RIGHT)

    statusbar.pack(side=BOTTOM, fill=X)

    # Log frame #
    logFrame = Frame(root, height=100, bd=0, highlightthickness=0, bg='white')
    logFrame.pack(side=BOTTOM, fill=X, expand=False)
    logFrame.pack_propagate(0)

    scrollbar = Scrollbar(logFrame, bg='#393c43', highlightbackground='#393c43', troughcolor='#393c43')
    scrollbar.pack(side=RIGHT, fill=Y)

    listbox = Listbox(logFrame, fg='#a9b7c6', bg='#393c43', bd=0, selectbackground='#262626', highlightthickness=0,
                      yscrollcommand=scrollbar.set)
    listbox.pack(side=LEFT, fill=BOTH, expand=True)
    scrollbar.config(command=listbox.yview)

    FTIR_fittingtool_GUI(root, root, listbox, statusbar, status1, status2)

    listbox.delete(0, END)
    listbox.insert(END, '*' * 60)
    listbox.insert(END, 'Welcome to FTIR Fitting Tool!')
    listbox.insert(END, 'This is the log file.')
    listbox.insert(END, 'Click to copy to the clipboard.')
    listbox.insert(END, '*' * 60)
    listbox.yview(END)

    root.mainloop()


if __name__ == '__main__':
    main()
