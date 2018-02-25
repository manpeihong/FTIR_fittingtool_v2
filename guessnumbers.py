import math
from tkinter import *
from random import randint
from sys import platform as _platform
import configparser
from ColorTheme import color_theme

__version__ = '0.4'


class guessnumbers_GUI(Frame):
    def __init__(self, root, masterroot, listbox):
        super().__init__(root, width=1000, bg='#2b2b2b')
        self.root = root
        self.masterroot = masterroot
        self.listbox = listbox
        self.digit = 4
        self.numberoftries = 10
        self.numberoftriesleft = self.numberoftries
        self.number = 0
        self.numberstring = ''
        self.myguess = ''
        self.As = 0
        self.Bs = 0
        self.exclude = []

        self.config = configparser.ConfigParser()
        self.config.read('configuration.ini')
        self.config_theme = self.config["Settings"]["colortheme"]
        self.bg = color_theme(self.config_theme).bg
        self.fg = color_theme(self.config_theme).fg
        self.bg_toolbar = color_theme(self.config_theme).bg_toolbar
        self.facecolor = color_theme(self.config_theme).facecolor
        self.warningcolor1 = color_theme(self.config_theme).warningcolor1
        self.warningcolor2 = color_theme(self.config_theme).warningcolor2
        self.warningcolor3 = color_theme(self.config_theme).warningcolor3

        self.frame0 = Frame(self, width=1000, bg=self.bg)
        self.frame0.pack(side=TOP, fill=X, expand=True)
        self.frame0.pack_propagate(0)

        buttonstart = Button(self.frame0, text="(S)tart",
                            command=self.startgame, highlightbackground=self.bg, width=10)
        buttonstart.grid(row=0, column=0, columnspan=1)

        Label(self.frame0, text='# of digits:',
              fg=self.fg, bg=self.bg, width=8, anchor=E).grid(row=0, column=1, sticky=E)
        self.entry_1 = Entry(self.frame0, highlightbackground=self.bg, width=8)
        self.entry_1.grid(row=0, column=2)
        self.entry_1.insert(0, '4')

        Label(self.frame0, text='# of tries:',
              fg=self.fg, bg=self.bg, width=8, anchor=E).grid(row=0, column=3, sticky=E)
        self.entry_2 = Entry(self.frame0, highlightbackground=self.bg, width=8)
        self.entry_2.grid(row=0, column=4)
        self.entry_2.insert(0, '10')

        Label(self.frame0, text='Guess a number:',
              fg=self.fg, bg=self.bg, width=13, anchor=E).grid(row=0, column=5, sticky=E)
        self.entry_3 = Entry(self.frame0, highlightbackground=self.bg, width=8)
        self.entry_3.grid(row=0, column=6)
        self.entry_3.insert(0, '')

        buttonenter = Button(self.frame0, text="Enter",
                            command=self.enternumber, highlightbackground=self.bg, width=10)
        buttonenter.grid(row=0, column=7, columnspan=1)

        if _platform == "darwin":
            buttonenter.config(text="Enter(⏎)")

        self.trylabel = Label(self.frame0, text='Tries left:', fg=self.fg, bg=self.bg, width=43, anchor=W)
        self.trylabel.grid(row=0, column=8, sticky=W)

        def buttonstart_event(event):
            self.startgame()

        def enternumber_event(event):
            self.enternumber()

        if _platform == "darwin":
            masterroot.bind('<Command-s>', buttonstart_event)
            masterroot.bind('<Return>', enternumber_event)
        elif _platform == "win32" or _platform == "win64" or _platform == "linux" or _platform == "linux2":
            masterroot.bind('<Control-s>', buttonstart_event)
            masterroot.bind('<Return>', enternumber_event)

        self.pack()

    def startgame(self):
        self.digit = int(self.entry_1.get())
        self.numberoftries = int(self.entry_2.get())
        self.numberoftriesleft = self.numberoftries
        self.randomnumber()
        # self.addlog(self.numberstring)
        self.addlog('Game begins!', self.warningcolor2)
        self.trylabel.config(text='Tries left: {}'.format(self.numberoftriesleft))

    def enternumber(self):
        self.myguess = str(self.entry_3.get())
        self.As = 0
        self.Bs = 0
        self.exclude = []
        if len(self.numberstring) == 0:
            self.startgame()

        if len(self.myguess) != self.digit:
            self.addlog('Invalid input.')
            return
        for i in range(0, self.digit):
            if self.myguess[i] == self.numberstring[i]:
                self.As += 1
        if self.As == self.digit:
            self.numberoftriesleft -= 1
            self.trylabel.config(text='Tries left: {}'.format(self.numberoftriesleft))
            self.addlog('Bingo! The number is {}.'.format(self.numberstring), self.warningcolor1)
            return
        for i in range(0, self.digit):
            for j in range(0, self.digit):
                if j not in self.exclude:
                    if self.myguess[j] == self.numberstring[i]:
                        self.Bs += 1
                        self.exclude.append(j)
                        break
        self.Bs -= self.As
        self.addlog('{}    {}A{}B'.format(self.myguess, self.As, self.Bs))
        self.numberoftriesleft -= 1
        self.trylabel.config(text='Tries left: {}'.format(self.numberoftriesleft))
        if self.numberoftriesleft <= 0:
            self.addlog('Game over! The number is {}'.format(self.numberstring), self.warningcolor2)

        self.entry_3.delete(0, END)

    def randomnumber(self):
        while True:
            oknumber = 1
            self.number = randint(0, math.pow(10, self.digit)-1)
            self.numberstring = "%0{}d".format(self.digit) % self.number
            for i in range(0, self.digit):
                for j in range(0, self.digit):
                    if j != i:
                        if self.numberstring[j] == self.numberstring[i]:
                            oknumber = 0
            if oknumber == 1:
                break

    def addlog(self, string, fgcolor="default"):
        self.listbox.insert(END, string)
        self.listbox.yview(END)

        if fgcolor != "default":
            i = 0
            while True:
                try:
                    self.listbox.itemconfig(i, bg=color_theme(self.config_theme).bg_log)
                    i += 1
                except:
                    self.listbox.itemconfig(i - 1, fg=fgcolor)
                    break


def main():
    root = Tk()
    w = 1000  # width for the Tk root
    h = 150  # height for the Tk root
    ws = root.winfo_screenwidth()  # width of the screen
    hs = root.winfo_screenheight()  # height of the screen
    x = (ws / 2) - (w / 2)
    y = (hs / 4) - (h / 4)
    root.geometry('%dx%d+%d+%d' % (w, h, x, y))

    root.wm_title("Guess A Number! v. {}".format(__version__))
    root.configure(background='#2b2b2b')

    # Status bar #
    statusbar = Frame(root, bg='#2b2b2b', bd=1, relief=RIDGE)

    authorLabel = Label(statusbar, text='© 10,2017 Peihong Man.', fg='#a9b7c6', bg='#2b2b2b', bd=1, relief=RIDGE,
                        padx=4.2, width=21)
    authorLabel.pack(side=LEFT)
    authorLabel.pack_propagate(0)
    status1 = Label(statusbar, text='Guess A Number!', fg='#a9b7c6', bg='#2b2b2b', bd=1, relief=RIDGE)
    status1.pack(side=LEFT, fill=X, expand=True)
    status1.pack_propagate(0)
    status2 = Label(statusbar, text='v. {}'.format(__version__), fg='#a9b7c6', bg='#2b2b2b', bd=1, relief=RIDGE, width=21)
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

    guessnumbers_GUI(root, root, listbox)

    listbox.delete(0, END)
    listbox.insert(END, '*' * 60)
    listbox.insert(END, 'Welcome to Guess A Number!')
    listbox.insert(END, '*' * 60)
    listbox.yview(END)

    root.mainloop()


if __name__ == '__main__':
    main()
