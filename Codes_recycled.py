def clearalldata_old(self):
    """Clear everything. """

    if self.layernumber != 0 or self.filepath.cget('text') != '':
        self.FTIRplot.clear()
        try:
            self.absorptionplot.clear()
        except AttributeError:
            pass
        self.lowercut = 400
        self.highercut = 6000
        self.transcut = 70
        self.transcutlow = 0
        self.FTIRplot.set_xlim([self.lowercut, self.highercut])
        self.FTIRplot.set_ylim([self.transcutlow, self.transcut])
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
        self.entry_21.insert(0, "0.7")
        self.entry_22.insert(0, "0.0")
        self.entry_23.insert(0, "0.0")
        self.entry_24.insert(0, "0.0")

        self.entry_31.delete(0, END)
        self.entry_32.delete(0, END)
        self.entry_33.delete(0, END)
        self.entry_332.delete(0, END)
        self.entry_31.insert(0, self.lowercut)
        self.entry_32.insert(0, self.highercut)
        self.entry_33.insert(0, self.transcutlow)
        self.entry_332.insert(0, self.transcut)

        self.layernumber = 0
        self.Temp = 300
        self.displayreflection, self.displayabsorption = 0, 0
        self.layertype_list, self.entry_x_list, self.entry_d_list, self.checklayer_list = [], [], [], []

        self.frame3.pack_forget()
        self.frame3 = Frame(self.frame3_shell, width=150, bg='#2b2b2b')
        self.frame3.pack(side=RIGHT, fill=BOTH, expand=True)
        self.frame3.pack_propagate(1)

        self.layertypevar1, self.layertypevar2, self.layertypevar3, self.layertypevar4, self.layertypevar5, \
        self.layertypevar6, self.layertypevar7, self.layertypevar8, self.layertypevar9, self.layertypevar10, \
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

        def change_sub(*args):
            if self.varsub.get() == "CdTe/Si":
                self.subtype = 1
            elif self.varsub.get() == "CdZnTe":
                self.subtype = 2

        Label(self.frame3, text='Layer:', fg="#a9b7c6", bg='#2b2b2b', width=self.COLUMN0_WIDTH).grid(row=10,
                                                                                                     column=0,
                                                                                                     sticky=E)
        Label(self.frame3, text='x:', fg="#a9b7c6", bg='#2b2b2b', width=self.COLUMN1_WIDTH + 1).grid(row=10,
                                                                                                     column=1,
                                                                                                     sticky=E)
        Label(self.frame3, text='d:', fg="#a9b7c6", bg='#2b2b2b', width=self.COLUMN2_WIDTH).grid(row=10, column=2,
                                                                                                 sticky=E)

        self.varsub = StringVar(self.frame3)
        self.varsub.set("Si")  # initial value
        self.varsub.trace("w", change_sub)
        suboption1 = OptionMenu(self.frame3, self.varsub, "Si", "CdZnTe")
        suboption1.config(bg="#2b2b2b", width=self.COLUMN0_WIDTH, anchor=E)
        suboption1.grid(row=11, column=0, sticky=W + E)

        self.entry_d_0 = Entry(self.frame3, highlightbackground='#2b2b2b', width=self.COLUMN2_WIDTH)
        self.entry_d_0.grid(row=11, column=2)
        self.entry_d_0.insert(0, '500')

        Label(self.frame3, text='  ', fg="#a9b7c6", bg='#2b2b2b', width=self.COLUMN0_WIDTH).grid(row=11, column=3,
                                                                                                 sticky=E)

        self.numberofdata = 0
        self.numberofdata2 = 0
        self.addlog('*' * 60)

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

def cal_trans(self):
    self.R23 = self.R2 + self.R3 * (1 - self.R2) * (1 - self.R2) * self.a2 * self.a2 / (
    1 - self.R2 * self.R3 * self.a2 * self.a2)
    self.T23 = (1 - self.R2) * (1 - self.R3) * self.a2 / (1 - self.R2 * self.R3 * self.a2 * self.a2)

    self.T13 = (1 - self.R1) * (1 - self.H) * self.T23 * self.a1 \
               / (1 - self.R1 * (1 - self.H) * self.R23 * self.a1 * self.a1)
    self.addlog(self.T13)
    return self.T13

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
