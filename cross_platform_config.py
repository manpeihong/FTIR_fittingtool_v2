from sys import platform as _platform

__version__ = '0.01'


class config:
    if _platform == "darwin" or _platform == "linux" or _platform == "linux2":
        # MAC OS X or linux
        FRAME_WIDTH = 1000
        FRAME_HEIGHT = 650
        LABEL_WIDTH = 10
        ENTRY_WIDTH = 8
        OPTIONMENU_WIDTH = 38

    elif _platform == "win32" or _platform == "win64":
        # Windows
        FRAME_WIDTH = 1150
        FRAME_HEIGHT = 920
        LABEL_WIDTH = 12
        ENTRY_WIDTH = 8
        OPTIONMENU_WIDTH = 34
