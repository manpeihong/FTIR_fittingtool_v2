

class color_theme:
    def __init__(self, theme):
        self.theme = int(theme)
        self.bg = ''
        self.fg = ''
        self.bg_toolbar = ''
        self.bg_log = ''
        self.facecolor = ''
        self.warningcolor1 = ''
        self.warningcolor2 = ''
        self.warningcolor3 = ''

        self.changetheme()

    def changetheme(self):
        if self.theme == 0:  # Dark: Default theme
            self.bg = "#2b2b2b"
            self.fg = "#a9b7c6"
            self.bg_toolbar = '#262626'
            self.bg_log = '#393c43'
            self.facecolor = 'gainsboro'
            self.warningcolor1 = 'red'
            self.warningcolor2 = 'yellow'
            self.warningcolor3 = 'royalblue'
        elif self.theme == 1:  # Light
            self.bg = "whitesmoke"
            self.fg = "dimgrey"
            self.bg_toolbar = 'gainsboro'
            self.bg_log = 'darkgrey'
            self.facecolor = 'white'
            self.warningcolor1 = 'red'
            self.warningcolor2 = 'yellow'
            self.warningcolor3 = 'royalblue'
        elif self.theme == 2:  # Royal
            self.bg = "orangered"
            self.fg = "white"
            self.bg_toolbar = 'gold'
            self.bg_log = 'gold'
            self.facecolor = 'lemonchiffon'
            self.warningcolor1 = 'red'
            self.warningcolor2 = 'green'
            self.warningcolor3 = 'royalblue'

        elif self.theme == 3:  # Sky
            self.bg = "lightskyblue"
            self.fg = "white"
            self.bg_toolbar = 'dodgerblue'
            self.bg_log = 'cornflowerblue'
            self.facecolor = 'lightcyan'
            self.warningcolor1 = 'red'
            self.warningcolor2 = 'yellow'
            self.warningcolor3 = 'darkviolet'

        elif self.theme == 4:  # Creamy
            self.bg = "papayawhip"
            self.fg = "gray"
            self.bg_toolbar = 'orange'
            self.bg_log = 'khaki'
            self.facecolor = 'seashell'
            self.warningcolor1 = 'red'
            self.warningcolor2 = 'green'
            self.warningcolor3 = 'royalblue'

        elif self.theme == 5:  # Mystery
            self.bg = "mediumpurple"
            self.fg = "white"
            self.bg_toolbar = 'darkviolet'
            self.bg_log = 'mediumorchid'
            self.facecolor = 'lavender'
            self.warningcolor1 = 'lime'
            self.warningcolor2 = 'orange'
            self.warningcolor3 = 'deepskyblue'

        elif self.theme == 6:  # Spring
            self.bg = "greenyellow"
            self.fg = "black"
            self.bg_toolbar = 'limegreen'
            self.bg_log = 'mediumseagreen'
            self.facecolor = 'palegreen'
            self.warningcolor1 = 'red'
            self.warningcolor2 = 'orange'
            self.warningcolor3 = 'royalblue'