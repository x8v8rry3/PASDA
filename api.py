# -*- coding: utf-8 -*-
from __future__  import division
from Tkinter import Tk
import os
from interface import *
		
root = Tk()
a = StartFace(root)
root.resizable(0,0)
root.mainloop()
root = Tk()
BasicDesk(root,int(a.value))
root.resizable(0,0)
root.mainloop()
