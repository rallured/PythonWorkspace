import matplotlib
matplotlib.use('TkAgg')

from numpy import arange, sin, pi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

import Tkinter as Tk
import sys

root = Tk.Tk()
root.wm_title("Embedding in TK")


f = Figure(figsize=(5,4), dpi=100)
a = f.add_subplot(211)
b = f.add_subplot(212)
t = arange(0.0,3.0,0.01)
s = sin(2*pi*t)

a.plot(t,s)


# a tk.DrawingArea
canvas = FigureCanvasTkAgg(f, master=root)
canvas.show()
canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

#Frame for Labels
frm = Tk.Frame(master=root)
frm.pack(side=Tk.BOTTOM,fill=Tk.BOTH,expand=1)

#toolbar = NavigationToolbar2TkAgg( canvas, root )
#toolbar.update()
canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

def _quit():
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate

def _test():
    b.plot(2*sin(2*pi*t))
    canvas.show()

button = Tk.Button(master=root, text='Quit', command=_quit)
button.pack(side=Tk.BOTTOM)
button2 = Tk.Button(master=root,text='Times2', command=_test)
button2.pack(side=Tk.BOTTOM)

Tk.mainloop()
