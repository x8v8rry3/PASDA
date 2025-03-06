# -*- coding: utf-8 -*-
from __future__  import division
from Tkinter import *
import ttk
from PIL import Image, ImageTk
import tkFileDialog
import tkMessageBox
import copy
from classes import *
from memory_pic import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import base64
import os

PeachHotcake = Sampling()
CatSandwich = Operation()
func0 = [u'([])ME',u'([])PBIAS',u'(-)MAE',u'(-)MSE',u'(-)RMSE',u'(-)RSR',u'(+)NSE',u'(+)R2',u'(+)R-S',u'(+)KGE']
value0 = [('me','[]'),('pbias','[]'),('mae','-'),('mse','-'),('rmse','-'),('rsr','-'),('ns','+'),('r2','+'),('r-s','+'),('kge','+')]
W = 600
H = 800

class StartFace(object):
	def __init__(self,main):
		self.mainw = main
		self.mainw.geometry('600x400')
		self.mainw.title('choose process')
		self.vartype = StringVar()
		self.vartype.set('1')
		Radiobutton(self.mainw, text='sampling', variable=self.vartype, font=('SimHei', 15), value='1').pack(pady=10)	
		Radiobutton(self.mainw, text='analysis', variable=self.vartype, font=('SimHei', 15), value='3').pack(pady=10)
		Radiobutton(self.mainw, text='sampling & analysis', variable=self.vartype, font=('SimHei', 15), value='4').pack(pady=10)
		Button(self.mainw, text='confirm', font=('SimHei', 15), command=self.live).pack(pady=10)

	def live(self):
		self.value = self.vartype.get()
		self.mainw.destroy()

class BasicDesk(object):
	def __init__(self,main,startfrom = 4,w=W,h=H/2):
		self.mainw = main
		self.mainw.title(u'PASDA(PArameter Sampling and Data Analysis toolkit)')
		self.mainw.geometry('%dx%d'%(w,h))
		self.mainw.startfrom = startfrom
		if startfrom == 1 or startfrom == 4:		
			FirstFace(self.mainw)	
		elif startfrom == 3:
			ThirdFace(self.mainw)
		self.__build_file_tree()
		
	def __build_file_tree(self):
		self.__make('./INPUT/SAMP_IN')
		self.__make('./OUTPUT/SAMP_OUT')
		self.__make('./INPUT/CALI_IN/OBS')
		self.__make('./INPUT/CALI_IN/SIM')
		self.__make('./INPUT/CALI_IN/PAR')
		self.__make('./OUTPUT/CALI_OUT/TEXT')
		self.__make('./OUTPUT/CALI_OUT/PIC')
		self.__make('./SAMPLE_DATA/INPUT/SAMP_IN')
		self.__make('./SAMPLE_DATA/INPUT/CALI_IN/OBS')
		self.__make('./SAMPLE_DATA/INPUT/CALI_IN/SIM')
		self.__make('./SAMPLE_DATA/INPUT/CALI_IN/PAR')
		if not os.path.exists('./INPUT/SAMP_IN/input.csv'):
			ff = open('./INPUT/SAMP_IN/input.csv','w')
			ff.close()
			
		get_pic(sim_csv,'./SAMPLE_DATA/INPUT/CALI_IN/SIM/sta1+FLOW.csv')
		get_pic(obs_csv,'./SAMPLE_DATA/INPUT/CALI_IN/OBS/sta1+FLOW.csv')
		get_pic(par_input,'./SAMPLE_DATA/INPUT/CALI_IN/PAR/par_input.txt')
		get_pic(sample_input,'./SAMPLE_DATA/INPUT/SAMP_IN/input.csv')		
		
				
	def __make(self,path):
		if not os.path.exists(path):
			os.makedirs(path)		
		

class FirstFace(object):
	def __init__(self,main):
		self.mainw = main
		self.firstface = Frame(height=H/2, width=W)
		self.firstface.pack_propagate(0) 
		self.firstface.pack()
		
		l = Label(self.firstface, text='get input data', bg='white', font=('SimHei', 25), width=30, height=2)
		l.pack()
		l = Button(self.firstface, text='>>', font=('SimHei', 12), width=5, height=2, command=self.next)
		l.place(relx=0.75, rely=0.7, anchor='s')
		b = Button(self.firstface, text='open input CSV file', font=('SimHei', 12), width=30, height=2, command=self.__open_csv)
		b.place(relx=0.35, rely=0.45, anchor='s')		
		b = Button(self.firstface, text='import data from CSV', font=('SimHei', 12), width=30, height=2, command=self.__do_read_csv)
		b.place(relx=0.35, rely=0.7, anchor='s')
		b = Button(self.firstface, text='view imported data', font=('SimHei', 12), width=30, height=2, command=self.__show_info)
		b.place(relx=0.35, rely=0.95, anchor='s')
	
	def next(self):
		isDone = PeachHotcake.show_par_info(doInfo = 0)
		if not isDone:
			stillContinue = tkMessageBox.askyesno(title=u'question', message=u'Parameter information has not been read yet, do you want to continue?', parent=self.firstface)
		else:
			stillContinue = 1
		if stillContinue:
			self.firstface.pack_forget()
			SecondFace(self.mainw)

	def __open_csv(self):
		os.startfile(PeachHotcake.inputpath)
		
	def __do_read_csv(self):
		string = PeachHotcake.read_input_csv(PeachHotcake.inputpath)
		if string == 'done':
			tkMessageBox.showinfo(title=u'information', message=u'read finish', parent=self.firstface)
		else: 
			tkMessageBox.showerror(title=u'error!', message=string, parent=self.firstface)
	
	def __show_info(self):
		string = PeachHotcake.show_par_info()
		wsii = Toplevel()
		wsii.resizable(0,0)
		wsii.attributes("-toolwindow", 1)
		wsii.attributes('-topmost',1)
		wsii.title(u'parameter detail')
		ls = string.split('\n')
		makeScroll(wsii,ls)

class SecondFace(object):	
	def __init__ (self,main):
		self.mainw = main
		self.secondface = Frame(height=H/2, width=W)
		self.secondface.pack_propagate(0)
		self.secondface.pack()
		
		l = Label(self.secondface, text='choose sampling way', bg='white', font=('SimHei', 25), width=30, height=2)
		l.pack()		
		b = Button(self.secondface, text='Random', font=('SimHei', 12), width=20, height=2, command=lambda:self.__do_sampling_main('simple'))
		b.place(relx=0.2, rely=0.54, anchor='s')
		b = Button(self.secondface, text='Latin-Hypercube', font=('SimHei', 12), width=20, height=2, command=lambda:self.__do_sampling_main('LHS'))
		b.place(relx=0.2, rely=0.86, anchor='s')
		b = Button(self.secondface, text='Quasi-Monte Carlo', font=('SimHei', 12), width=20, height=2, command=self.__quasiAsk)
		b.place(relx=0.6, rely=0.54, anchor='s')
		b = Button(self.secondface, text='MCMC', font=('SimHei', 12), width=20, height=2, command=self.__mcmcAsk)
		b.place(relx=0.6, rely=0.86, anchor='s')

		if self.mainw.startfrom == 1:
			b = Button(self.secondface, text='<<', font=('SimHei', 12), width=5, height=2, command=self.back)
			b.place(relx=0.9, rely=0.7, anchor='s')
		
		else:
			b = Button(self.secondface, text='>>', font=('SimHei', 12), width=5, height=2, command=self.next)
			b.place(relx=0.9, rely=0.54, anchor='s')
			b = Button(self.secondface, text='<<', font=('SimHei', 12), width=5, height=2, command=self.back)
			b.place(relx=0.9, rely=0.86, anchor='s')
		
	def next(self):
		self.secondface.pack_forget()
		ThirdFace(self.mainw)		
		
	def back(self):
		self.secondface.pack_forget()
		FirstFace(self.mainw)

	def __quasiAsk(self):
		wsii = Toplevel()
		wsii.resizable(0,0)
		wsii.attributes("-toolwindow", 1)
		wsii.attributes('-topmost',1)
		wsii.title(u'Quasi Monte Carlo sampling parameters')
		wsii.geometry('300x300')
		vartype = StringVar()
		vartype.set('halton')
		varplace = IntVar()
		
		l1 = Label(wsii, text=u'sampling type:', font=('SimHei', 12), width=30)
		r1 = Radiobutton(wsii, text='Halton', variable=vartype, font=('SimHei', 15), value='halton')
		r2 = Radiobutton(wsii, text='Faure', variable=vartype, font=('SimHei', 15), value='faure')
		r3 = Radiobutton(wsii, text='Sobol', variable=vartype, font=('SimHei', 15), value='sobol')
		l2 = Label(wsii, text=u'Sequence truncation position', font=('SimHei', 12), width=30)
		e1 = Entry(wsii, show=None, font=('SimHei', 15),textvariable=varplace)
		b1 = Button(wsii, text=u'conform', font=('SimHei', 15), width=30, command = lambda: self.__GoOn(vartype,varplace,wsii,'quasi') )
		l1.pack()
		r1.pack()
		r2.pack()
		r3.pack()
		l2.pack()
		e1.pack()
		b1.pack()

	def __mcmcAsk(self):
		wsii = Toplevel()
		wsii.resizable(0,0)
		wsii.attributes("-toolwindow", 1)
		wsii.attributes('-topmost',1)
		wsii.title(u'MCMC sampling parameters')
		wsii.geometry('300x200')
		vartype = StringVar()
		vartype.set('mcmc')
		
		l1 = Label(wsii, text=u'sampling type:', font=('SimHei', 12), width=30)
		r1 = Radiobutton(wsii, text='MCMC', variable=vartype, font=('SimHei', 15), value='mcmc')
		r2 = Radiobutton(wsii, text='M-H', variable=vartype, font=('SimHei', 15), value='mh')
		r3 = Radiobutton(wsii, text='Gibbs', variable=vartype, font=('SimHei', 15), value='gibbs')
		b1 = Button(wsii, text=u'conform', font=('SimHei', 15), width=30, command = lambda: self.__GoOn(vartype,0,wsii,'mcmc') )
		l1.pack()
		r1.pack()
		r2.pack()
		r3.pack()
		b1.pack()
		
		
	def __GoOn(self,p1,p2,window,way):
		try:
			if way == 'quasi':
				PeachHotcake.quasitype = p1.get()
				PeachHotcake.cutdown = p2.get()				
				self.__do_sampling_main('quasi')
				window.destroy()
			elif way == 'mcmc':
				PeachHotcake.mcmctype = p1.get()				
				self.__do_sampling_main('MCMC')
				window.destroy()
		except ValueError:
			tkMessageBox.showerror(title=u'error!', message=u"Sequence truncation position is not integer please reinput", parent=window)
		except TypeError:
			tkMessageBox.showerror(title=u'error!', message=u"Without uniform distribution parameters, low difference sequence sampling cannot be performed", parent=window)
			
			
	def __do_sampling_main(self,way): 
		output, namels= PeachHotcake.do_sampling(way)
		if namels == ['quasiwrong!']:
			raise TypeError
		ls = list(output)
		outputname = 'output_' + way + '.txt'
		outstr = 'NAME'
		tempup = 'UP'
		tempdown = 'DOWN'
		string = ''
		for ii in namels:
			outstr = outstr + '\t' + PeachHotcake.parls[ii].name
			tempup = tempup + '\t' + str(PeachHotcake.parls[ii].up)
			tempdown = tempdown + '\t' + str(PeachHotcake.parls[ii].down)
		if way == 'quasi':
			string = u'WARNING! This sampling method only supports uniform distribution, and other distribution parameters will be skipped\n'
			outstr = outstr + '\t' + PeachHotcake.quasitype + '\n'
		elif way == 'MCMC':
			outstr = outstr + '\t' + PeachHotcake.mcmctype + '\n'
		else:
			outstr = outstr + '\t' + way + '\n'
		outstr = outstr + tempup + '\n' + tempdown + '\n'
		count = 0
		for sample in ls:
			count += 1
			ls0 = []
			for value in sample:
				ls0.append(str(value))
			outstr = outstr + str(count) + '\t' + '\t'.join(ls0) + '\n'
			string = string + u'sampling paraseter set%d：'% count
			if namels == []:
				namels = range(len(ls0))
			for ii in range(len(namels)): 
				string = string + u'parameter %s = %f，'%(PeachHotcake.parls[namels[ii]].name,float(ls0[ii]))
			string = string + '\n'
		string = string + u'sampling results are saved in SAMP_OUT/' + unicode(outputname) +u'中\n-------------------\n'

		with open('./OUTPUT/SAMP_OUT/' + outputname,'w') as ff:
			ff.write(str(outstr))	
		wsii = Toplevel()
		wsii.resizable(0,0)
		wsii.attributes("-toolwindow", 1)
		wsii.attributes('-topmost',1)
		wsii.title(u'sampling detail')
		ls = string.split('\n')
		makeScroll(wsii,ls,font=('SimHei', 14))
		
class ThirdFace(object):	
	def __init__ (self,main):
		self.mainw = main
		self.thirdface = Frame(height=H/2, width=W)
		self.thirdface.pack_propagate(0)
		self.thirdface.pack()		
		
		l = Label(self.thirdface, text='model analysis-input', bg='white', font=('SimHei', 25), width=30, height=2)
		l.pack()			
		b = Button(self.thirdface, text='read parameter', font=('SimHei', 12), width=25, height=2, command=self.__read_par)
		b.place(relx=0.35, rely=0.45, anchor='s')		
		b = Button(self.thirdface, text='read simulation and observation', font=('SimHei', 12), width=40, height=2, command=self.__read_sim_obs)
		b.place(relx=0.35, rely=0.7, anchor='s')
		b = Button(self.thirdface, text='read caculation result', font=('SimHei', 12), width=25, height=2, command=self.__read_calc)
		b.place(relx=0.35, rely=0.95, anchor='s')	
		if self.mainw.startfrom == 3:
			b = Button(self.thirdface, text='>>', font=('SimHei', 12), width=5, height=2, command=self.next)
			b.place(relx=0.75, rely=0.7, anchor='s')
		else:
			b = Button(self.thirdface, text='<<', font=('SimHei', 12), width=5, height=2, command=self.back)
			b.place(relx=0.75, rely=0.86, anchor='s')	
			b = Button(self.thirdface, text='>>', font=('SimHei', 12), width=5, height=2, command=self.next)
			b.place(relx=0.75, rely=0.54, anchor='s')

	def back(self):
		self.thirdface.pack_forget()
		SecondFace(self.mainw)
	
	def next(self):
		self.thirdface.pack_forget()
		FourthFace(self.mainw)	
	
	def __read_par(self): 
		CatSandwich.erase()
		string = CatSandwich.read_par('./INPUT/CALI_IN/PAR/par_input.txt')
		if string == 'done':
			tkMessageBox.showinfo(title=u'information', message=u'FINISH', parent=self.thirdface)
			CatSandwich.nosim = False
		else: 
			tkMessageBox.showerror(title=u'error!', message=string, parent=self.thirdface)
		
	def __read_sim_obs(self): 
		if CatSandwich.nosim:
			tkMessageBox.showinfo(title=u'error', message=u'please read parameter', parent=self.thirdface)
			return 0
		ls = os.listdir('./INPUT/CALI_IN/SIM')		
		allperfect = 1
		for name in ls:
			string = CatSandwich.read_sim_obs_csv('./INPUT/CALI_IN/SIM/' + name, 'sim')
			if string != 'done':
				allperfect = 0
				tkMessageBox.showerror(title=u'error!', message=string, parent=self.thirdface)
		if allperfect:
			tkMessageBox.showinfo(title=u'information', message=u'simulation read done', parent=self.thirdface)			
		ls = os.listdir('./INPUT/CALI_IN/OBS')
		allperfect = 1
		for name in ls:
			string = CatSandwich.read_sim_obs_csv('./INPUT/CALI_IN/OBS/' + name, 'obs')
			if string != 'done':
				allperfect = 0
				tkMessageBox.showerror(title=u'error!', message=string, parent=self.thirdface)
		if allperfect:
			tkMessageBox.showinfo(title=u'information', message=u'observation read done', parent=self.thirdface)	  
			
	def __read_calc(self):
		path = tkFileDialog.askopenfilename(title=u'choose file')
		string = CatSandwich.read_file(path)
		if string == 'done':
			tkMessageBox.showinfo(title=u'information', message=u'parameter read finish', parent=self.thirdface)
		else:
			tkMessageBox.showerror(title=u'error!', message=string, parent=self.thirdface)			
			
class FourthFace(object):	
	def __init__ (self,main):
		self.mainw = main	
		self.fourthface = Frame(height=H/2, width=W)
		self.fourthface.pack_propagate(0)
		self.fourthface.pack()	
		l = Label(self.fourthface, text='model analysis set', bg='white', font=('SimHei', 25), width=30, height=2)
		l.pack()			
		b = Button(self.fourthface, text='select sensitivity analysis method', font=('SimHei', 12), width=45, height=2, command=self.__sensiAsk)
		b.place(relx=0.35, rely=0.45, anchor='s')		
		b = Button(self.fourthface, text='select uncertainty analysis method', font=('SimHei', 12), width=45, height=2, command=self.__uncerAsk)
		b.place(relx=0.35, rely=0.7, anchor='s')
		b = Button(self.fourthface, text='selsect performance criteria', font=('SimHei', 12), width=45, height=2, command=self.__likeAsk)
		b.place(relx=0.35, rely=0.95, anchor='s')		
		b = Button(self.fourthface, text='<<', font=('SimHei', 12), width=5, height=2, command=self.back)
		b.place(relx=0.75, rely=0.86, anchor='s')	
		b = Button(self.fourthface, text='>>', font=('SimHei', 12), width=5, height=2, command=self.next)
		b.place(relx=0.75, rely=0.54, anchor='s')	

	def __sensiAsk(self): 
		wsii = Toplevel()
		wsii.resizable(0,0)
		wsii.attributes("-toolwindow", 1)
		wsii.attributes('-topmost',1)
		wsii.title(u'sensitivity method')
		wsii.geometry('300x200')
		vartype = StringVar()
		vartype.set('mlr')
		
		l1 = Label(wsii, text=u'method', font=('SimHei', 12), width=30)
		r1 = Radiobutton(wsii, text=u'MLR', variable=vartype, font=('SimHei', 15), value='mlr')
		r2 = Radiobutton(wsii, text=u'S-K test', variable=vartype, font=('SimHei', 15), value='k-s')
		r3 = Radiobutton(wsii, text=u'Sobol test', variable=vartype, font=('SimHei', 15), value='sobol')
		r4 = Radiobutton(wsii, text=u'skip', variable=vartype, font=('SimHei', 15), value='0')
		b1 = Button(wsii, text=u'conform', font=('SimHei', 15), width=30, command = lambda: self.__GoOn(vartype,wsii,'sensi') )
		l1.pack()
		r1.pack()
		r2.pack()
		r3.pack()
		r4.pack()
		b1.pack()

	def __uncerAsk(self):
		wsii = Toplevel()
		wsii.resizable(0,0)
		wsii.attributes("-toolwindow", 1)
		wsii.attributes('-topmost',1)
		wsii.title(u'uncertainty method')
		wsii.geometry('300x200')
		vartype = StringVar()
		vartype.set('hsy')
		
		l1 = Label(wsii, text=u'method', font=('SimHei', 12), width=30)
		r1 = Radiobutton(wsii, text=u'HSY', variable=vartype, font=('SimHei', 15), value='hsy')
		r2 = Radiobutton(wsii, text=u'GLUE', variable=vartype, font=('SimHei', 15), value='glue')
		r3 = Radiobutton(wsii, text=u'SUFI-2', variable=vartype, font=('SimHei', 15), value='sufi2')
		r4 = Radiobutton(wsii, text=u'skip', variable=vartype, font=('SimHei', 15), value='0')
		b1 = Button(wsii, text=u'conform', font=('SimHei', 15), width=30, command = lambda: self.__GoOn(vartype,wsii,'uncer'))
		l1.pack()
		r1.pack()
		r2.pack()
		r3.pack()
		r4.pack()
		b1.pack()		
	
	def __likeAsk(self): 
		wsii = Toplevel()
		wsii.resizable(0,0)
		wsii.attributes("-toolwindow", 1)
		wsii.attributes('-topmost',1)
		wsii.title(u' performance criteria and accept threshold')
		wsii.geometry('600x400')
		var = []
		threshold = []
		
		frame_l = Frame(wsii,height=400,width=300)
		frame_r = Frame(wsii,height=400,width=300)
		frame_b = Frame(wsii,height=50,width=600)
		frame_l.pack_propagate(0)
		frame_r.pack_propagate(0)
		frame_b.pack_propagate(0)
		frame_l.place(relx=0,rely=0,anchor='nw')
		frame_r.place(relx=1,rely=0,anchor='ne')
		frame_b.place(relx=0,rely=1,anchor='sw')
		
		self.func = copy.deepcopy(func0)
		self.value = copy.deepcopy(value0)
		if CatSandwich.likedic and CatSandwich.simlist[0].stalist[0].sim == None:
			for tup in zip(self.func,self.value):
				if tup[1][0] not in list(CatSandwich.likedic):
					self.func.remove(tup[0])
					self.value.remove(tup[1])
		if CatSandwich.uncertype == 'glue':
			self.func = self.func[6:10]
			self.value = self.value[6:10]
	
		l1 = Label(frame_l, text=u'performance criteria', font=('SimHei', 12), width=30)
		l1.place(relx=-0.05,rely=0,anchor='nw')
		ii = 0
		for name in self.func:
			intVar = IntVar()
			intVar.set(0)
			var.append(intVar)
			Checkbutton(frame_l, text=name, variable=intVar, onvalue=ii+1, font=('SimHei', 15)).place(relx=0.15,rely=0.1+0.07*ii)
			ii += 1
		
		l2 = Label(frame_r, text=u'accept threshold', font=('SimHei', 12), width=30)
		l2.place(relx=-0.1,rely=0,anchor='nw')
		jj = 0
		for name in self.func:
			if name == u'([])ME' or name == u'([])PBIAS':
				doubleVar1 = DoubleVar()
				doubleVar2 = DoubleVar()
				threshold.append([doubleVar1,doubleVar2])
				Entry(frame_r, show=None, font=('SimHei', 15),textvariable=doubleVar1).place(relx=0.15,rely=0.11+0.07*jj,relwidth=0.37)
				Entry(frame_r, show=None, font=('SimHei', 15),textvariable=doubleVar2).place(relx=0.53,rely=0.11+0.07*jj,relwidth=0.37)
			else:
				doubleVar = DoubleVar()
				threshold.append(doubleVar)
				Entry(frame_r, show=None, font=('SimHei', 15),textvariable=doubleVar).place(relx=0.15,rely=0.11+0.07*jj,relwidth=0.75)
			jj += 1
		
		b1 = Button(frame_b, text=u'conform', font=('SimHei', 15), width=30, command = lambda: self.__GoOn([var,threshold],wsii,'like') )
		b1.pack(side = 'top')
		
	def __GoOn(self,p1,window,way):
		if way == 'sensi':	
			window.destroy()
			CatSandwich.sensitype = p1.get()	
		elif way == 'uncer':		
			window.destroy()
			CatSandwich.uncertype = p1.get()		
		elif way == 'like':
			try:
				funcls = [v.get()-1 for v in p1[0]]
				thls = []	
				for th in p1[1]:
					if isinstance(th, list):
						temp = [th[0].get(),th[1].get()]
						thls.append(temp)
					else:
						thls.append(th.get())
				dic = {}							
				tkMessageBox.showinfo(title=u'information', message=u'done', parent=window)				
				for (func,th) in zip(funcls, thls):
					if func >= 0:
						func = self.value[func]
						dic[func[0]] = [th,func[1]]
				CatSandwich.likedic = dic
				window.destroy()
			except ValueError:
				tkMessageBox.showerror(title=u'error!', message=u"threshold wrong reinput", parent=window)
	
	def back(self):
		self.fourthface.pack_forget()
		ThirdFace(self.mainw)
	
	def next(self):
		self.fourthface.pack_forget()
		FifthFace(self.mainw)				

		
class FifthFace(object):	
	def __init__ (self,main):
		self.mainw = main
		self.fifthface = Frame(height=H/2, width=W)
		self.fifthface.pack_propagate(0)
		self.fifthface.pack()		
		
		l = Label(self.fifthface, text='analysis-calculation & perfoem', bg='white', font=('SimHei', 15), width=45, height=2)
		l.pack()			
		b = Button(self.fifthface, text='calculate perfomation criteria', font=('SimHei', 12), width=45, height=2, command=self.calculation)
		b.place(relx=0.35, rely=0.45, anchor='s')
		b = Button(self.fifthface, text='sensitivity and uncertainty analysis', font=('SimHei', 12), width=45, height=2, command=self.askpostin)
		b.place(relx=0.35, rely=0.7, anchor='s')
		b = Button(self.fifthface, text='show result', font=('SimHei', 12), width=45, height=2, command= self.show)
		b.place(relx=0.35, rely=0.95, anchor='s')	
		b = Button(self.fifthface, text='<<', font=('SimHei', 12), width=5, height=2, command=self.back)
		b.place(relx=0.75, rely=0.7, anchor='s')	

	def calculation(self):
		string = CatSandwich.do_calculation()
		if string == 'done':
			tkMessageBox.showinfo(title=u'information', message=u'finish')
		CatSandwich.write_output()
	
	def askpostin(self):
		wsii = Toplevel()
		wsii.resizable(0,0)
		wsii.attributes("-toolwindow", 1)
		wsii.attributes('-topmost',1)		
		wsii.title(u'choose station')
		wsii.geometry('600x400')
		
		frame_1 = Frame(wsii,height=300,width=100)
		frame_2 = Frame(wsii,height=300,width=100)
		frame_3 = Frame(wsii,height=300,width=100)
		frame_b = Frame(wsii,height=50,width=100)
		frame_1.pack_propagate(0)
		frame_2.pack_propagate(0)
		frame_3.pack_propagate(0)
		frame_b.pack_propagate(0)
		frame_1.place(relx=0.2,rely=0,anchor='nw')
		frame_2.place(relx=0.5,rely=0,anchor='n')
		frame_3.place(relx=0.8,rely=0,anchor='ne')		
		frame_b.place(relx=0,rely=1,anchor='sw')
		
		stalist = []
		var = StringVar()
		like = StringVar()
		var.set(0)
		like.set(0)
		Label(frame_1, text=u'station', font=('SimHei', 15), width=50).pack()
		Label(frame_2, text=u'variable', font=('SimHei', 15), width=50).pack()
		Label(frame_3, text=u'method', font=('SimHei', 15), width=50).pack()
		for name in CatSandwich.stanamels:
			stringvar = StringVar()
			stringvar.set(0)
			stalist.append(stringvar)
			Checkbutton(frame_1, text=name, variable=stringvar, onvalue=name, font=('SimHei', 15), width=10).pack()
		for name in CatSandwich.varnamels:
			Radiobutton(frame_2, text=name, variable=var, font=('SimHei', 15), value=name, width=10).pack()
		for name in list(CatSandwich.likedic):
			if CatSandwich.uncertype == 'glue' and name not in ['ns','r2','r-s','kge']:
				pass
			else:
				Radiobutton(frame_3, text=name, variable=like, font=('SimHei', 15), value=name, width=10).pack()
		b1 = Button(frame_b, text=u'conform', font=('SimHei', 15), width=10, command = lambda: self.__Do(wsii,[[v.get() for v in stalist if v.get() != '0'],var.get(),like.get()])).pack()
		

	def __Do(self,window,postin):		
		window.destroy()
		CatSandwich.postin = postin
		try:
			CatSandwich.do_uncertainty()
			CatSandwich.do_sensitivity()
			tkMessageBox.showinfo(title=u'information', message=u'FINISH')
		except:
			tkMessageBox.showerror(title=u'error!', message=u'FAILURE')
		try:
			CatSandwich.do_display()		
			tkMessageBox.showinfo(title=u'information', message=u'FINISH')
		except:
			tkMessageBox.showerror(title=u'error!', message=u'FAILURE')
		
	def show(self):
		wsii = Toplevel()
		wsii.resizable(0,0)
		wsii.attributes('-toolwindow', 1)
		wsii.attributes('-topmost',1)		
		wsii.title(u'analysis result')
		wsii.geometry('1000x600')
		
		tabControl = ttk.Notebook(wsii)	
		tab1 = Frame(tabControl)
		tabControl.add(tab1,text=u'sensitivity result')
		tab2 = Frame(tabControl)
		tabControl.add(tab2,text=u'uncertainty-parameter')
		if CatSandwich.simlist[0].stalist[0].obs:
			tab3 = Frame(tabControl)
			tabControl.add(tab3,text=u'uncertainty-variable')
		
		tab1_f1 = Frame(tab1, height=630, width=400)
		tab1_f2 = Frame(tab1, height=630, width=600)
		tab1_f3 = Frame(tab1, height=70, width=1000)
		tab1_f1.pack_propagate(0)
		tab1_f2.pack_propagate(0)
		tab1_f3.pack_propagate(0)
		tab1_f1.place(relx=0, rely=0, anchor='nw')
		tab1_f2.place(relx=1, rely=0, anchor='ne')
		tab1_f3.place(relx=1, rely=1, anchor='se')
		columns = (u'parname',u'p-value',u'p-value')
		Label(tab1_f1, text='result-sheet', font=('SimHei', 15), width=35).pack(side=TOP,pady=10)
		style = ttk.Style(tab1_f1)
		style.configure('Treeview', rowheight=15)
		tab1_tv = ttk.Treeview(tab1_f1, height=30, show='headings', columns=columns)
		tab1_tv.column(u'parname', width=120, anchor='center')
		tab1_tv.column(u'p-value', width=90, anchor='center')
		tab1_tv.column(u'p-value', width=90, anchor='center')
		tab1_tv.heading(u'parname', text=u'parname')
		tab1_tv.heading(u'p-value', text=u'p-value')
		tab1_tv.heading(u'p-value', text=u'p-value')
		tab1_tv.pack(side=TOP,pady=5)		
		for i in range(len(CatSandwich.s_result)):
			tab1_tv.insert('',i,values=(CatSandwich.parnamels[i],CatSandwich.s_result[i][0],CatSandwich.s_result[i][1]))		

		Label(tab1_f2, text='result-bar', font=('SimHei', 15), width=25).pack(side=TOP,pady=10)
		self.__show_pic(tab1_f2,CatSandwich.sensi_pic, side=TOP, fill=BOTH, pady=5)
		b1 = Button(tab1_f3, text=u'save sheet & figure', font=('SimHei', 15), width=20, command = lambda:self.__save_sensi(wsii)).pack()		
		tab2_f1 = Frame(tab2, height=630, width=480)
		tab2_f2 = Frame(tab2, height=630, width=480)
		tab2_f3 = Frame(tab2, height=70, width=1000)
		tab2_f1.pack_propagate(0)
		tab2_f2.pack_propagate(0)
		tab2_f3.pack_propagate(0)
		tab2_f1.place(relx=0.01, rely=0, anchor='nw')
		tab2_f2.place(relx=0.99, rely=0, anchor='ne')
		tab2_f3.place(relx=1, rely=1, anchor='se') 
		Label(tab2_f1, text=u'parametervalue  scatter', font=('SimHei', 15), width=25).pack(side=TOP,pady=10)			
		Button(tab2_f1, text=u'next parameter', font=('SimHei', 15), width=20, command = lambda:self.__show_pic_handler(tab2_f1,1,side=TOP, fill=BOTH, pady=0)).pack(side=TOP, pady=5)
		self.tab2_c1 = self.__show_pic(tab2_f1,CatSandwich.par_plot_picls[0],side=TOP, fill=BOTH, pady=0)
		Label(tab2_f2, text=u'parvalue distribution histogram', font=('SimHei', 15), width=25).pack(side=TOP,pady=10)		
		Button(tab2_f2, text=u'next parameter', font=('SimHei', 15), width=20, command = lambda:self.__show_pic_handler(tab2_f2,2,side=TOP, fill=BOTH, pady=0)).pack(side=TOP, pady=5)		
		self.tab2_c2 = self.__show_pic(tab2_f2,CatSandwich.posterior_picls[0],side=TOP, fill=BOTH, pady=0)
		Button(tab2_f3, text=u'save sheet and figure', font=('SimHei', 15), width=20, command = lambda:self.__save_uncerpar(wsii)).place(x=500,y=0,anchor='n')	
		if CatSandwich.simlist[0].stalist[0].obs:
			tab3_f1 = Frame(tab3, height=325, width=680)
			tab3_f2 = Frame(tab3, height=325, width=680)
			tab3_f3 = Frame(tab3, height=325, width=320)
			tab3_f4 = Frame(tab3, height=325, width=320)
			tab3_f5 = Frame(tab3, height=50, width=1000)
			tab3_f1.pack_propagate(0)
			tab3_f2.pack_propagate(0)
			tab3_f3.pack_propagate(0)
			tab3_f4.pack_propagate(0)
			tab3_f5.pack_propagate(0)
			tab3_f1.place(relx=0, rely=0, anchor='nw')
			tab3_f2.place(relx=0, rely=325/700, anchor='nw')
			tab3_f3.place(relx=1, rely=0, anchor='ne')	
			tab3_f4.place(relx=1, rely=325/700, anchor='ne')
			tab3_f5.place(relx=1, rely=1, anchor='se')	
			Label(tab3_f1, text='variable time series figure', font=('SimHei', 15), width=25).pack(side=TOP,pady=10)
			Button(tab3_f1, text=u'next station', font=('SimHei', 15), width=20, command = lambda:self.__show_pic_handler(tab3_f1,3,side=TOP, fill=BOTH, pady=0)).pack(side=TOP, pady=5)		
			self.tab3_c1 = self.__show_pic(tab3_f1,CatSandwich.time_pic[0],side=TOP, fill=BOTH, pady=0)			
			Label(tab3_f2, text='best sim', font=('SimHei', 15), width=25).pack(side=TOP,pady=10)
			Button(tab3_f2, text=u'next station', font=('SimHei', 15), width=20, command = lambda:self.__show_pic_handler(tab3_f2,4,side=TOP, fill=BOTH, pady=0)).pack(side=TOP, pady=5)	
			self.tab3_c2 = self.__show_pic(tab3_f2,CatSandwich.time_best[0],side=TOP, fill=BOTH, pady=0)
			Label(tab3_f3, text='uncertainty factor', font=('SimHei', 15), width=25).pack(side=TOP,pady=10)
			columns = (u'station name',u'p-factor',u'r-factor')
			tab3_tv1 = ttk.Treeview(tab3_f3, height=15, show='headings', columns=columns)
			tab3_tv1.column(u'station name', width=100, anchor='center')
			tab3_tv1.column(u'p-factor', width=90, anchor='center')
			tab3_tv1.column(u'r-factor', width=90, anchor='center')
			tab3_tv1.heading(u'station name', text=u'station name')
			tab3_tv1.heading(u'p-factor', text=u'p-factor')
			tab3_tv1.heading(u'r-factor', text=u'r-factor')
			tab3_tv1.pack(side=TOP,pady=5)		
			for i in range(len(CatSandwich.p_factor)):
				tab3_tv1.insert('',i,values=(CatSandwich.postin[0][i],CatSandwich.p_factor[i],CatSandwich.r_factor[i]))				
			Label(tab3_f4, text='best parameter set', font=('SimHei', 15), width=25).pack(side=TOP,pady=10)
			columns = (u'parname1',u'parvalue1',u'parname2',u'parvalue2')
			tab3_tv2 = ttk.Treeview(tab3_f4, height=15, show='headings', columns=columns)
			tab3_tv2.column(u'parname1', width=75, anchor='center')
			tab3_tv2.column(u'parvalue1', width=75, anchor='center')
			tab3_tv2.column(u'parname2', width=75, anchor='center')
			tab3_tv2.column(u'parvalue2', width=75, anchor='center')
			tab3_tv2.heading(u'parname1', text=u'parname')
			tab3_tv2.heading(u'parvalue1', text=u'parvalue')
			tab3_tv2.heading(u'parname2', text=u'parname')
			tab3_tv2.heading(u'parvalue2', text=u'parvalue')
			tab3_tv2.pack(side=TOP,pady=5)		
			for i in range(0,len(CatSandwich.parnamels),2):
				if i+1 < len(CatSandwich.parnamels)-1:
					tab3_tv2.insert('',i,values=(CatSandwich.parnamels[i],CatSandwich.bestpar[i],CatSandwich.parnamels[i+1],CatSandwich.bestpar[i+1]))
				else:
					tab3_tv2.insert('',i,values=(CatSandwich.parnamels[i],CatSandwich.bestpar[i],'',''))	
			Button(tab3_f5, text=u'save sheet & figure', font=('SimHei', 15), width=20, command = lambda:self.__save_uncervar(wsii)).place(x=500,y=0,anchor='n')
		tabControl.pack(expand=1, fill="both")		
		
	def __show_pic(self,tab,pic,**kw):	
		canvas = FigureCanvasTkAgg(pic, master=tab)
		canvas.draw()
		canvas.get_tk_widget().pack(**kw)
		return canvas
		
	def __show_pic_handler(self,tab,loc,**kw):
		if loc == 1:
			a = CatSandwich.par_plot_picls.pop(0)
			CatSandwich.par_plot_picls.append(a)
			pic = CatSandwich.par_plot_picls[0]
			self.tab2_c1.get_tk_widget().destroy()
			self.tab2_c1 = self.__show_pic(tab,pic=pic,**kw)
		elif loc == 2:
			a = CatSandwich.posterior_picls.pop(0)
			CatSandwich.posterior_picls.append(a)
			pic = CatSandwich.posterior_picls[0]
			self.tab2_c2.get_tk_widget().destroy()
			self.tab2_c2 = self.__show_pic(tab,pic=pic,**kw)
		elif loc == 3:
			a = CatSandwich.time_pic.pop(0)
			CatSandwich.time_pic.append(a)
			pic = CatSandwich.time_pic[0]
			self.tab3_c1.get_tk_widget().destroy()
			self.tab3_c1 = self.__show_pic(tab,pic=pic,**kw)		
		elif loc == 4:
			a = CatSandwich.time_best.pop(0)
			CatSandwich.time_best.append(a)
			pic = CatSandwich.time_best[0]
			self.tab3_c2.get_tk_widget().destroy()
			self.tab3_c2 = self.__show_pic(tab,pic=pic,**kw)				
				
	def __save_sensi(self,window):
		a = CatSandwich.save_pic(CatSandwich.sensi_pic,'s_result.png')
		b = CatSandwich.save_sensi()
		window.attributes('-topmost',0)
		if a == 'done' and b == 'done':
			tkMessageBox.showinfo(title=u'information', message=u'allsaved')
		elif a == 'not' and b == 'done':
			tkMessageBox.showinfo(title=u'WARNING', message=u'only file saved')
		elif a == 'done' and b == 'not':
			tkMessageBox.showinfo(title=u'WARNING', message=u'only figure saved')
		elif a == 'not' and b == 'not':
			tkMessageBox.showinfo(title=u'WARNING', message=u'save failure')
		window.attributes('-topmost',1)

	def __save_uncerpar(self,window):
		for ii in range(len(CatSandwich.parnamels)): 
			CatSandwich.save_pic(CatSandwich.par_plot_picls[ii], CatSandwich.parnamels[ii]+'_plot.jpg')
			CatSandwich.save_pic(CatSandwich.posterior_picls[ii], CatSandwich.parnamels[ii]+'_posterior.jpg')
		CatSandwich.save_parpost()
		window.attributes('-topmost',0)
		tkMessageBox.showinfo(title=u'information', message=u'saved')
		window.attributes('-topmost',1)
		
	def __save_uncervar(self,window):
		for ii in range(len(CatSandwich.postin[0])): 
			CatSandwich.save_pic(CatSandwich.time_pic[ii], CatSandwich.postin[0][ii]+'_95PPU.jpg')
			CatSandwich.save_pic(CatSandwich.time_best[ii], CatSandwich.postin[0][ii]+'_best.jpg')
		CatSandwich.save_varpost()
		window.attributes('-topmost',0)
		tkMessageBox.showinfo(title=u'information', message=u'saved')
		window.attributes('-topmost',1)
		
	def back(self):
		self.fifthface.pack_forget()
		FourthFace(self.mainw)
		
def makeScroll(obj,msg,font=('SimHei', 12),width=50,height=20):
	scrollbar1 = Scrollbar(obj)
	scrollbar2 = Scrollbar(obj,orient=HORIZONTAL)
	scrollbar1.pack(side='right', fill=Y)	
	scrollbar2.pack(side='bottom', fill=X, anchor=S)
	listbox = Listbox(obj, yscrollcommand=scrollbar1.set, xscrollcommand=scrollbar2.set, font=font, width=width, height=height)
	for line in msg:
		listbox.insert(END, line)
	listbox.pack(fill=BOTH, expand=YES)
	scrollbar1.config(command=listbox.yview)
	scrollbar2.config(command=listbox.xview)

def get_pic(pic_code, pic_name):
	image = open(pic_name, 'wb')
	image.write(base64.b64decode(pic_code))
	image.close()
