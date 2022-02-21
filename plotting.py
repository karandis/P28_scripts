#plotting.py
import matplotlib.pyplot as pl
#------------------------------- DEFAULT ---------------------------------------
# pl.style.use('seaborn-dark')
#Font and linewidth
pl.rcParams['font.weight'] = "bold"
pl.rcParams['font.family'] = "sans-serif"
pl.rcParams['font.size'] = 20
pl.rcParams['legend.fontsize'] = 20
pl.rcParams['lines.linewidth'] = 2
pl.rcParams['axes.linewidth'] = 3
pl.rcParams['figure.figsize'] = (12,8)
#Axes, label and face colors
pl.rcParams['axes.facecolor'] = 'white'
pl.rcParams['figure.facecolor'] = 'white'
pl.rcParams['patch.facecolor'] = 'white'
pl.rcParams['savefig.facecolor'] = 'white'
# pl.rcParams['savefig.bbox'] = 'tight'
pl.rcParams['axes.labelcolor'] = 'black'
pl.rcParams['axes.edgecolor'] = 'black'
pl.rcParams['figure.edgecolor'] = 'black'
pl.rcParams['ytick.color']= 'black'	 # color
pl.rcParams['xtick.color']= 'black'	 # color
#axes
pl.rcParams['axes.grid'] = False
pl.rcParams['axes.titlesize'] = 15
pl.rcParams['axes.labelsize'] = 18
pl.rcParams['axes.labelweight'] = 'medium'
pl.rcParams['axes.spines.left'] = True
pl.rcParams['axes.spines.right'] = True
pl.rcParams['axes.spines.top'] = True
pl.rcParams['axes.spines.bottom'] = True
#ticks
pl.rcParams['xtick.top']= False	 # xtick on top
pl.rcParams['xtick.bottom']= True	 # xtick on top
pl.rcParams['xtick.minor.visible']= True
pl.rcParams['xtick.labelsize']= 15	 # xtick labelsize size in points
pl.rcParams['xtick.major.size']= 12	 # major tick size in points
pl.rcParams['xtick.minor.size']= 8	 # minor tick size in points
pl.rcParams['xtick.direction']= 'out'	 # in, out or inout
pl.rcParams['ytick.left']= True	 # ytick left
pl.rcParams['ytick.minor.visible']= True	 # ytick left
pl.rcParams['ytick.right']= False	 # ytick right
pl.rcParams['ytick.labelsize']= 15	 # xtick labelsize size in points
pl.rcParams['ytick.major.size']= 12	 # major tick size in points
pl.rcParams['ytick.minor.size']= 8	 # minor tick size in points
pl.rcParams['ytick.direction']= 'out'	 # in, out or inout
# --------------------------- BEGIN USER INPUT ---------------------------------

def set_colors(facecolor,edgecolor):
    pl.rcParams['axes.facecolor'] = facecolor
    pl.rcParams['figure.facecolor'] = facecolor
    pl.rcParams['legend.facecolor'] = facecolor
    pl.rcParams['patch.facecolor'] = facecolor
    pl.rcParams['savefig.facecolor'] = facecolor
    pl.rcParams['axes.labelcolor'] = edgecolor
    pl.rcParams['axes.edgecolor'] = edgecolor
    pl.rcParams['legend.edgecolor'] = edgecolor
    pl.rcParams['text.color'] = edgecolor
    pl.rcParams['figure.edgecolor'] = edgecolor
    pl.rcParams['ytick.color']= edgecolor
    pl.rcParams['xtick.color']= edgecolor

def set_mode(x):
    if x == 'paper':
        #---------------------------- For papers -------------------------------
        #Font and linewidth
        pl.rcParams['font.weight'] = "medium"
        pl.rcParams['font.family'] = "sans-serif"
        pl.rcParams['font.size'] = 17
        pl.rcParams['legend.fontsize'] = 14
        pl.rcParams['lines.linewidth'] = 1.5
        pl.rcParams['axes.linewidth'] = 1.5
        pl.rcParams['figure.figsize'] = (8,6)
        #Axes, label and face colors
        pl.rcParams['axes.facecolor'] = 'white'
        pl.rcParams['figure.facecolor'] = 'white'
        pl.rcParams['patch.facecolor'] = 'white'
        pl.rcParams['savefig.facecolor'] = 'white'
        # pl.rcParams['savefig.bbox'] = 'tight'
        pl.rcParams['axes.labelcolor'] = 'black'
        pl.rcParams['axes.edgecolor'] = 'black'
        pl.rcParams['figure.edgecolor'] = 'black'
        pl.rcParams['ytick.color']= 'black'	 # color
        pl.rcParams['xtick.color']= 'black'	 # color
        #axes
        pl.rcParams['axes.grid'] = False
        pl.rcParams['axes.titlesize'] = 18
        pl.rcParams['axes.labelsize'] = 17
        pl.rcParams['axes.labelweight'] = 'medium'
        pl.rcParams['axes.spines.left'] = True
        pl.rcParams['axes.spines.right'] = True
        pl.rcParams['axes.spines.top'] = True
        pl.rcParams['axes.spines.bottom'] = True
        #ticks
        pl.rcParams['xtick.top']= False	 # xtick on top
        pl.rcParams['xtick.bottom']= True	 # xtick on top
        pl.rcParams['xtick.minor.visible']= True
        pl.rcParams['xtick.labelsize']= 12	 # xtick labelsize size in points
        pl.rcParams['xtick.major.size']= 10	 # major tick size in points
        pl.rcParams['xtick.minor.size']= 7	 # minor tick size in points
        pl.rcParams['xtick.direction']= 'out'	 # in, out or inout
        pl.rcParams['ytick.left']= True	 # ytick left
        pl.rcParams['ytick.minor.visible']= True	 # ytick left
        pl.rcParams['ytick.right']= False	 # ytick right
        pl.rcParams['ytick.labelsize']= 12	 # xtick labelsize size in points
        pl.rcParams['ytick.major.size']= 10	 # major tick size in points
        pl.rcParams['ytick.minor.size']= 7	 # minor tick size in points
        pl.rcParams['ytick.direction']= 'out'	 # in, out or inout
