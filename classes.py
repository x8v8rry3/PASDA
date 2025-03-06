# -*- coding: utf-8 -*-
from __future__  import division
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import csv
import sobol_seq 
import copy
import itertools
from common import *
from quasibasic import dim_max, dim_vector
from scipy.special import comb
from scipy.stats import ks_2samp, pearsonr, spearmanr, kendalltau, rankdata, t
 
class Parameter(object): 
	def __init__(self,name,cdf):
		self.name = name
		self.cdf = cdf
		self.pdf = self.build_pdf(self.cdf)
		self.down = cdf[0][0]
		self.up = cdf[-1][0]
		
	def build_pdf(self,cdf): 
		cdf.sort(key = lambda x: x[0])
		pdf = {} 
		for ii in range(len(cdf)-1):
			key = (cdf[ii][0], cdf[ii+1][0]) 
			value = round(cdf[ii+1][1] - cdf[ii][1],6)
			pdf[key] = value
		return pdf
		
	def generate_info(self):
		str0 = u'参数名称：' + self.name + '\n'
		ls = list(self.pdf)
		ls.sort(key = lambda x: x[0])
		str1 = u'参数上界：' + str(ls[-1][1]) + '\n'
		str2 = u'参数下界：' + str(ls[0][0]) +  '\n'
		str3 = u'参数分布：\n'
		for key in ls:
			if key[1] == ls[-1][1]:
				str3 = str3 + u'    P(%.3f≤%s≤%.3f) = %f'%(key[0],self.name,key[1],self.pdf[key])
			else:
				str3 = str3 + u'    P(%.3f≤%s＜%.3f) = %f'%(key[0],self.name,key[1],self.pdf[key])
			str3 = str3 + '\n'
		string = str0 + str1 + str2 + str3 + '-----------------------'
		return string		
			
	def linear_mapping(self, direction = 'obverse',ls = []): 
		if direction == 'obverse':
			ls0 = self.pdf.keys()
			ls0.sort()
			ls = ls * (ls0[-1][1] - ls0[0][0]) + ls0[0][0]
			ls = np.round(ls,6)
			return ls
		elif direction == 'reverse':
			cdf = copy.deepcopy(self.cdf)
			delta = cdf[-1][0] - cdf[0][0]
			low = cdf[0][0]
			for ii in range(len(cdf)):
				cdf[ii] = (round((cdf[ii][0] - low)/delta,6),cdf[ii][1])			
			pdf = self.build_pdf(cdf)
			return pdf
			
class Sampling(object):

	def __init__(self, path = os.path.abspath('.')+'/INPUT/SAMP_IN/input.csv'):
		self.inputpath = path
		self.parls = []
		self.quasitype = 'sobol'
		self.mcmctype = 'gibbs'
		self.cutdown = 0
		
	def read_input_csv(self,fullpath):	
		with open(fullpath,'r') as ff:
			reader = csv.reader(ff)
			listreader = list(reader)
			if len(listreader[0]) != 5:
				return u"sample time input error：column number is not 5，Please refer to the examples and instructions to input parameters "			
			test = listreader[0][-1].decode('GB18030')
			if not test.isdecimal():
				return u"sample time input error：sample time not integer"
			if int(test) < 1:
				return u"sample time input error：sample time <1"	
			self.sampletime = int(test)
			
			parStart = False 
			cdf = []
			for ls in listreader:
				if ls[0]:				
					parStart = bool(1-parStart)
					if not parStart:
						if not self.strictly_increasing([x[0] for x in cdf]):
							return u"par%sCumulative probability distribution input error: Interval segmentation points are non strictly monotonically increasing"%name
						if not self.strictly_increasing([x[1] for x in cdf]):
							return u"par%sCumulative probability distribution input error: Cumulative probability non strictly monotonically increasing"%name
						if cdf[0][1] != 0 or cdf[-1][1] != 1:
							return u"par%sCumulative probability distribution input error: Cumulative probability upper and lower bounds error"%name													
						par = Parameter(name,cdf)
						cdf = []
						self.parls.append(par)
						parStart = bool(1-parStart)
					name = ls[0]
					if name == u'END':
						break
				if parStart:
					try:
						cdf.append((float(ls[1]),float(ls[2])))
					except ValueError:
						return u"Cumulative probability distribution input error: There are non floating point numbers present"
			return 'done'		
					
	def strictly_increasing(self,ls):	
		return all(x<y for x,y in zip(ls,ls[1:]))
		
	def show_par_info(self,doInfo = 1):
		strls = []
		if len(self.parls) == 0:
			string = u'\nno parameter,\nno data\n'
			isDone = False
		else:
			isDone = True
			string = ''
			if doInfo:
				for par in self.parls:
					strls.append(par.generate_info())
				string = '\n'.join(strls)
				string = u'sample time:' + str(self.sampletime) + '\n-*-*-*-*-*-*-*-*-*-*-*-\n' + string			
		if doInfo:
			return string
		else:
			return isDone
		
	def do_sampling(self,way):
		pdfls = []
		namels = range(len(self.parls))
		for par in self.parls:
			pdfls.append(par.linear_mapping(direction = 'reverse'))
		if way == 'LHS':
			result = self.latin_hypercube_sampling(pdfls)
		elif way == 'simple':
			result = self.simple_sampling(pdfls)		
		elif way == 'quasi':
			result, namels = self.quasi_sampling(pdfls)
		elif way == 'MCMC':
			result = self.mcmc_sampling(pdfls)
		if namels == ['quasiwrong!']:
			return ([], ['quasiwrong!'])
		resultarr = np.array(result)
		samplels = []
		count = 0
		for ii in namels:
			par = self.parls[ii]
			samplels.append(par.linear_mapping(direction = 'obverse', ls = resultarr[:,count]))
			count += 1
		samplels = np.array(samplels).T
		return samplels, namels

	def latin_hypercube_sampling(self, distri = []):
		parnum = len(self.parls)
		sampletime = self.sampletime
		result = np.empty([sampletime,parnum])
		dd = 1.0 / sampletime
		if distri == []:
			distri = [{(0,1):1}]*parnum
		save = []
		for dict in distri:
			distrils = sorted(dict.items(),key = lambda d:d[0][0])
			dvdpoint = self.__divide(distrils,sampletime)
			save.append(dvdpoint)
		for ii in range(len(save)):
			secls = self.__rebuildls(save[ii])
			savesample = []
			for section in secls:
				samp = np.random.uniform(section[0],section[1])
				savesample.append(samp)
			temp = np.array(savesample)
			np.random.shuffle(temp)
			result[:,ii] = temp.T
		return result
		
	def simple_sampling(self, distri = []):
		parnum = len(self.parls)
		sampletime = self.sampletime
		result = np.empty([sampletime,parnum])
		if distri == []:
			distri = [{(0,1):1}]*parnum
			
		for ii in range(parnum):
			distrils = sorted(distri[ii].items(),key = lambda d:d[0][0])
			section = [x[0] for x in distrils]
			pb = [x[1] for x in distrils]
			placesamp = np.random.choice(len(distrils),sampletime,p=pb) 
			temp = map(lambda x: np.random.uniform(section[x][0],section[x][1]), placesamp) 
			result[:,ii] = np.array(temp).T
		return result
		
	def quasi_sampling(self, distri = []):
		parnum = len(self.parls)
		sampletime = self.sampletime
		if distri == []:
			distri = [{(0,1):1}]*parnum
		namels = []
		for ii in range(len(distri)):
			if len(distri[ii]) == 1:
				namels.append(ii)
		distri = [x for x in distri if len(x) == 1]
		ndims = len(namels)
		if ndims == 0:
			return ([], ['quasiwrong!'])
		output = self.__sequence_sampling(ndims)
		return output, namels
		
	def mcmc_sampling(self, distri = {}):
		parnum = len(self.parls)
		sampletime = self.sampletime
		result = np.empty([sampletime,parnum])
		if distri == []:
			distri = [{(0,1):1}]*parnum
		namels = []
		
		if self.mcmctype == 'mcmc' or  self.mcmctype == 'mh': 
			for ii in range(parnum):
				distrils = sorted(distri[ii].items(),key = lambda d:d[0][0])
				section = [x[0] for x in distrils]
				pb = [x[1] for x in distrils]			
				placesamp =  self.__mcmc_mh(pb, sampletime) 
				temp = map(lambda x: np.random.uniform(section[x][0],section[x][1]), placesamp) 
				result[:,ii] = np.array(temp).T	
				
		elif self.mcmctype == 'gibbs': 
			allsec = []
			allpb = []
			for ii in range(parnum):
				distrils = sorted(distri[ii].items(),key = lambda d:d[0][0])
				section = [x[0] for x in distrils]
				pb = [x[1] for x in distrils]
				allsec.append(section)
				allpb.append(pb)
			placesamp =  self.__gibbs(allpb, parnum, sampletime) 
			arr = np.array(placesamp)
			temp = []
			for jj in range(parnum):
				temp.append(map(lambda x: np.random.uniform(allsec[jj][x][0],allsec[jj][x][1]), arr[:,jj]))
			result = np.array(temp).T
		return result
			
	def __mcmc_mh(self, p, sampletime, itertime=10000):
		stmat = np.array([p for x in range(len(p))], dtype=np.float64) 
		lsmc = [np.random.randint(len(p))] 
		lsoutput = []
		itercount = 0
		samplecount = 0
		while True:
			now = lsmc[samplecount] 
			next = np.random.choice(a=len(p),size=1,p=p)[0] 
			itercount += 1
			if self.mcmctype == 'mh':
				keep = (p[next]*stmat[next][now])/(p[now]*stmat[now][next]) 
				accept = min(keep,1)
			elif self.mcmctype == 'mcmc':
				keep = p[next]*stmat[next][now] 
			thold = np.random.uniform(0,1) 
			if thold < keep:
				samplecount += 1
				lsmc.append(next)
				if itercount > itertime:
					lsoutput.append(next)
			if len(lsoutput) >= sampletime:
				break
		return lsoutput
				
	def __gibbs(self, p, parnum, sampletime, itertime=10000):
		base = 0
		itercount = 0
		lsoutput = []
		lsmc = []
		first = []
		for kk in range(parnum):
			first.append(np.random.choice(a=len(p[kk]),size=1,p=p[kk])[0])
		lsmc.append(first)
		while True:
			next = lsmc[itercount] 
			new = copy.deepcopy(next)
			for ii in range(parnum):
				axis = base % parnum 
				if axis == 0:
					base = 0
				new[axis] = np.random.choice(a=len(p[axis]),size=1,p=p[axis])[0]
			itercount += 1
			lsmc.append(new)
			if itercount > itertime:
				lsoutput.append(new)			
			if len(lsoutput) >= sampletime:
				break				
			base += 1
		return lsoutput		
			
		
	def __divide(self,ls,secnum): 
		now = 0
		dd = round(1/secnum,6)
		ddbasic = dd
		lsback = [0]
		ddChanged = False
		for ii in range(len(ls)):
			sec = list(ls[ii][0])
			length = ls[ii][1]
			while dd < length:
				if now < secnum-1:
					point = self.__findpoint(sec,length,dd)
					now += 1
					lsback.append(point)
					length = length - dd
					sec[0] = point	
				else: 
					break
				if ddChanged:
					dd = ddbasic
					ddChanged = False
			dd = dd - length
			ddChanged = True
		lsback.append(1)
		return lsback
		
		
	def __findpoint(self,sec, length, aimlen): 
		point = round(sec[0]+(sec[1]-sec[0])*aimlen/length, 6)
		return point
		

	def __rebuildls(self,point): 
		ls = []
		for ii in range(len(point)-1):
			ls.append([point[ii],point[ii+1]])
		return ls

	def __sequence_sampling(self, ndims, random = 0):
		seqnum = self.sampletime
		cutdown = self.cutdown
		nsample = (ndims+2)*seqnum
		if random:
			i = int(np.floor(10*np.random.random()*seqnum))
		if cutdown:
			i = cutdown
		else:
			i = 10
		if self.quasitype == 'sobol':
			matrix = sobol_seq.i4_sobol_generate(2*ndims, seqnum+i)
		if self.quasitype == 'halton':
			base = dim_vector[0:2*ndims]
			matrix = np.zeros((seqnum+i,2*ndims))
			start = np.random.choice(range(100*ndims,300*ndims))
			ls = range(start, start+seqnum+i)	
			k = 0 					
			for num in ls:
				base_inv = 1 / np.array(base)
				t = np.ones(2*ndims)
				t = num * t
				ls = []
				while (0 < np.sum(t)):
					for j in range(0, 2*ndims):
						d = (t[j] % base[j])
						matrix[k,j] = matrix[k,j] + float(d) * base_inv[j]
						base_inv[j] = base_inv[j] / base[j]
						t[j] = (t[j] // base[j])
				k += 1	
				
		if self.quasitype == 'faure':
			j = 0
			while (dim_vector[j] < 2*ndims):
				j += 1
			base = dim_vector[j]
			matrix = np.zeros((seqnum+i,2*ndims))
			start = np.random.choice(range(100*ndims,300*ndims))
			k = 0
			while (k < seqnum+i):
				basenum = start
				base_inv = 1 / base
				ls = []
				while (0 < basenum):
					d = basenum % base
					ls.append(float(d))
					base_inv = base_inv / base
					basenum = basenum // base					
				matrix[k,0] = self.__calc_digital(ls, base)
				cmat = self.__build_transmatrix(len(ls))
				for jj in range(1,2*ndims):
					ls = np.matmul(cmat, np.array(ls).T)
					ls = ls % base
					matrix[k,jj] = self.__calc_digital(ls, base)
				k += 1
				start += 1		
		matrix = matrix[i:,:]
		matrix1 = matrix[:,0:ndims]
		matrix2 = matrix[:,ndims:]
		matrix = np.zeros(((ndims+2)*seqnum, ndims))
		for ii in range(0,ndims):
			matrixtemp = copy.deepcopy(matrix1)
			matrixtemp[:,ii] = matrix2[:,ii]
			matrix[(ii+2)*seqnum:(ii+3)*seqnum, :] = matrixtemp
		matrix[0:seqnum,:] = matrix1
		matrix[seqnum:2*seqnum,:] = matrix2
		return matrix
	
	def __calc_digital(self, ls, base):
		f = 0
		for ii in range(len(ls)):
			f += ls[ii] * ((1/base) ** (ii+1))
		return f
		
	def __build_transmatrix(self, ndim):
		mat = np.zeros((ndim, ndim))
		for ii in range(ndim):
			for jj in range(ii, ndim):
				mat[ii,jj] = comb(jj,ii)
		return mat

class Station(object):
	def __init__(self,name):
		self.staname = name
		self.obs = {}
		self.sim = {}
		self.calc = {}
		
	def back_calc_result(self,key = []):
		back = []
		for dicitem in self.calc.items():
			if key == []:
				back.extend(dicitem.values())
			else:
				temp = []
				for var in key:
					if var in dicitem:
						temp.append(dicitem[var])
				back.extend(temp)
		return back	
		
	def calc_do(self, way, var, sim, obs):
		if way == 'ns':
			self.calc[var]['ns'] = 1-np.sum((sim-obs)**2)/np.sum((obs-np.mean(obs))**2)
		elif way == 'pbias':
			self.calc[var]['pbias'] = np.sum(sim-obs)/np.sum(obs)
		elif way == 'r2':
			r,p = pearsonr(sim,obs)
			self.calc[var]['r2'] = r**2
		elif way == 'r-s': 
			if len(np.unique(sim)) == len(sim) and len(np.unique(obs)) == len(obs): 
				r,p = spearmanr(sim,obs)
				self.calc[var]['r-s'] = r
			else:
				a = rankdata(sim)
				b = rankdata(obs)
				r,p = pearsonr(a,b)
				self.calc[var]['r-s'] = r
		elif way == 'me':
			self.calc[var]['me'] = np.sum(sim-obs)/len(sim)
		elif way == 'mae':
			self.calc[var]['mae'] = np.sum(np.abs(sim-obs))/len(sim)
		elif way == 'mse':
			self.calc[var]['mse'] = np.sum((sim-obs)**2)/len(sim)
		elif way == 'rmse':
			self.calc[var]['rmse'] = (np.sum((sim-obs)**2)/len(sim))**0.5
		elif way == 'rsr':
			self.calc[var]['rsr'] = ((np.sum((sim-obs)**2))**0.5)/((np.sum((obs-np.mean(obs))**2))**0.5)
		elif way == 'kge':
			alpha = ((np.sum((sim-np.mean(sim))**2))**0.5)/((np.sum((obs-np.mean(obs))**2))**0.5)
			beta = np.mean(sim)/np.mean(obs)
			rele, p = pearsonr(sim,obs) 
			self.calc[var]['kge'] = 1-((rele-1)**2 + (alpha-1)**2 + (beta-1)**2)**0.5
		
class Simulation(object):
	def __init__(self,num):		
		self.pardict = {}
		self.stalist = []
		self.accept = 'yes'
		self.num = num
	
	def search_result(self,staname,varname,likename):
		if not isinstance(staname, list):
			staname = [staname]
		kk = 0
		for sta in self.stalist:
			if sta.staname in staname:
				kk += sta.calc[varname][likename]
		kk = kk/len(staname)
		return kk

	def do_calculation(self,key):
		ls = list(self.stalist[0].obs)
		for sta in self.stalist:
			if not sta.calc:
				sta.calc = dict.fromkeys(ls)
			for var in ls:
				sim = np.array(sta.sim[var])
				obs = np.array(sta.obs[var])
				if not sta.calc[var]:
					sta.calc[var] = {}
				sta.calc_do(key, var, sim, obs)

	def do_judgement(self,key,value):
		if self.accept == 'yes':
			for sta in self.stalist:
				for var in list(sta.calc):
					if value[1] == '+' and sta.calc[var][key] < value[0]:
						self.accept = 'no'
					elif value[1] == '-' and sta.calc[var][key] > value[0]:
						self.accept = 'no'
					elif value[1] == '[]' and sta.calc[var][key] > value[0][1]:
						self.accept = 'no'
					elif value[1] == '[]' and sta.calc[var][key] < value[0][0]:
						self.accept = 'no'
		
	def generate_result(self):				
		a = sorted(self.pardict.iteritems(), key=lambda d:d[0])
		a = [x[1] for x in a]
		for sta in self.stalist:
			for var in list(sta.calc):
				for func in list(sta.calc[var]):
					a.append(sta.calc[var][func])
		a.append(self.accept)
		return a
	
		
class Operation(object):
	def __init__(self):
		self.simlist = []
		self.simtime = 0
		self.nosim = True
		self.nosta = True
		self.sensetype = '0'
		self.uncertype = '0'
		self.likedic = {}
		self.stanamels = []
		self.varnamels = []
		self.parnamels = []
		self.parupdown = []
		self.postin = []
		
	def read_par(self,fullpath):
		self.erase(flag = True)
		with open(fullpath, 'r') as ff:
			lines = ff.readlines()
			self.parnamels = std_split(lines[0])[1:-1]
			tempup = std_split(lines[1])[1:]
			tempdown = std_split(lines[2])[1:]
			try:
				tempup = [float(x) for x in tempup]
				tempdown = [float(x) for x in tempdown]
			except ValueError:
				return u'存在不为浮点数的参数上下界值，请修改后重试'			
			self.parupdown = [[tempdown[ii],tempup[ii]] for ii in range(len(tempup))]
			for line in lines[3:]:
				try:
					ls = std_split(line, tofloat = True)
				except ValueError:
					return u'存在不为浮点数的参数值，请修改后重试'
				d = dict(zip(self.parnamels,ls[1:]))	
				if self.nosim:
					temp = Simulation(int(ls[0]))
					temp.pardict = d
					self.simlist.append(temp)
				else:
					self.simlist[int(ls[0])-1].pardict = d	
			self.parnamels.sort()
		return 'done'
		
	def read_sim_obs_csv(self,fullpath,way):
		reader = csv.reader(open(fullpath,'r'))
		name = os.path.split(fullpath)[1]
		staname = name.split('.')[0].split('+')[0]
		varname = name.split('.')[0].split('+')[1]
		if way == 'sim':
			if staname not in self.stanamels:
				self.nosta = True
				self.stanamels.append(staname)
			if varname not in self.varnamels:
				self.varnamels.append(varname)
		ii = 0
		for row in reader:
			try:
				ls = tuple([float(x) for x in row[1:]])
			except ValueError:
				return u'%s: 文件%s 行%d 存在不为浮点数的值'%(way, name, ii+1)
			if way == 'sim':
				if self.nosta:
					tempsta = Station(staname)
					tempsta.sim[str(varname)] = ls
					self.simlist[ii].stalist.append(tempsta)
				else:
					self.simlist[ii].stalist[-1].sim[str(varname)] = ls
			elif way == 'obs':
				for jj in range(len(self.simlist)):
					for kk in range(len(self.simlist[jj].stalist)):
						if self.simlist[jj].stalist[kk].staname == staname:
							self.simlist[jj].stalist[kk].obs[str(varname)] = ls	
			ii += 1
		if way == 'sim':
			self.simtime = ii
		self.nosta = False
		return 'done'
		
	def read_file(self,path):
		self.erase()
		with open(path,'r') as ff:
			reader = csv.reader(ff)
			ls = list(reader)
			for ii in range(len(ls)):
				ls[ii] = [x for x in ls[ii] if x]
			for ii in range(1,len(ls[5])):
				if '[' not in ls[6][ii]:
					try:
						temp = float(ls[6][ii])
					except ValueError:
						return u'行7单元%d存在不符合要求的值，请查看说明文档'%(ii+1)
					self.likedic[ls[5][ii]] = [temp,ls[7][ii]]
				else:
					try:
						a = float(ls[6][ii].split(',')[0][1:])
						b = float(ls[6][ii].split(',')[1][:-1])
					except ValueError:
						return u'行7单元%d存在不符合要求的值，请查看说明文档'%(ii+1)
					self.likedic[ls[5][ii]] = [[a,b],ls[7][ii]]
			self.parnamels = ls[0][1:]
			
			self.parupdown = []
			for ii in range(1,len(ls[0])):
				try:
					tempup = float(ls[1][ii])
				except ValueError:
					return u'行2单元%d存在非浮点数，请查看说明文档'%(ii+1)
				try:
					tempdown = float(ls[2][ii])
				except ValueError:
					return u'行3单元%d存在非浮点数，请查看说明文档'%(ii+1)
				self.parupdown.append([tempdown,tempup])
			self.stanamels = ls[3][1:]
			self.varnamels = ls[4][1:]
			self.simtime = 0
			lp = len(self.parnamels)
			ll = len(self.likedic)
			lst = len(self.stanamels)
			lv = len(self.varnamels)
			for jj in range(10,len(ls)):
				try:
					self.simtime += 1
					sim = Simulation(int(ls[jj][0]))
					par = [float(a) for a in ls[jj][1:(lp+1)]]
					sim.pardict = dict(zip(self.parnamels, par))
					sim.accept = ls[jj][-1]
					start = lp+1
					for kk in range(0,lst):
						sta = Station(self.stanamels[kk])
						sta.calc = dict.fromkeys(self.varnamels)
						for tt in range(0,lv):
							v = [float(a) for a in ls[jj][start:start+ll]]
							start += ll
							sta.calc[self.varnamels[tt]] = dict(zip(list(self.likedic),v))
						sim.stalist.append(sta)
					self.simlist.append(sim)
				except ValueError:
					return '模拟%d输入存在错误，请查看说明文档进行修改'%(jj-9)
		return 'done'

	def erase(self,flag = False):
		self.simlist = []
		self.simtime = 0
		self.nosim = True
		self.nosta = True
		self.likedic = {}
		self.stanamels = []
		self.varnamels = []
		self.parnamels = []
		self.parupdown = []
		self.postin = []
		if flag:
			self.sensitype = '0'
			self.uncertype = '0'

	def do_calculation(self):	
		for sim in self.simlist:
			sim.accept = 'yes'
			for key,value in self.likedic.items():
				sim.do_calculation(key)
				sim.do_judgement(key,value)
		return 'done'

	def do_uncertainty(self):
		if self.postin == []:
			self.postin = [self.stanamels[0],self.varnamels[0],list(self.likedic)[0]]
		self.simnum_a, pvalue_a, pvalue_r, like_a, best_par = self.accept_reject()
		self.likearr = np.array(like_a)
		self.pvalue_a = np.array(pvalue_a)
		self.pvalue_r = np.array(pvalue_r)
		self.bestpar = np.array(best_par)
		
		if self.uncertype == 'hsy':
			self.u_parpost = np.sort(self.pvalue_a, axis=0)	
			self.u_dispost = np.ones([self.pvalue_a.shape[0],len(self.parnamels)])/(self.pvalue_a.shape[0]+1)
			for ii in range(1,self.u_dispost.shape[0]):
				self.u_dispost[ii,:] = self.u_dispost[ii,:] + self.u_dispost[ii-1,:]
			
		elif self.uncertype == 'glue':
			postsort = np.argsort(self.pvalue_a, axis=0)
			self.u_parpost = np.sort(self.pvalue_a, axis=0)
			ppost = []
			for ii in range(len(self.parnamels)):
				ppost.append(self.likearr)
			ppost = np.array(ppost).T
			for jj in range(ppost.shape[1]):
				ppost[:,jj] = ppost[:,jj][postsort[:,jj]]
			self.u_dispost = (ppost/np.sum(ppost, axis=0))*(1-1/(ppost.shape[0]+1))
			for ii in range(1,self.u_dispost.shape[0]):
				self.u_dispost[ii,:] = self.u_dispost[ii,:] + self.u_dispost[ii-1,:]
			
		elif self.uncertype == 'sufi2':
			if not hasattr(self,'sufi2post'):
				self.u_parpost =  np.array(self.parupdown).T
			if self.simtime > 500:
				return u'采样次数超过500，请减少采样次数'
			j = np.zeros([int(comb(self.simtime,2)),len(self.parnamels)])
			count = 0
			pvalue = np.array(self.__generate_parlist())
			for ii in itertools.combinations(range(self.simtime),2):
				j[count,:] = (self.likearr[ii[0]]-self.likearr[ii[1]])/(pvalue[ii[0]]-pvalue[ii[1]])
				count += 1
			h = np.dot(j.T,j)
			c = np.std(self.likearr) ** 2 * np.linalg.inv(h)
			sj = np.diag(c) ** 0.5
			pdown = self.bestpar - sj * t.ppf(0.975, self.simtime-len(self.parnamels))
			pup = self.bestpar + sj * t.ppf(0.975, self.simtime-len(self.parnamels))
			x = (pdown-self.u_parpost[0,:])/2
			y = (self.u_parpost[1,:]-pup)/2
			tempdown = pdown - np.maximum(x,y)
			tempup = pup + np.maximum(x,y)
			self.u_parpost[0,:] = tempdown
			self.u_parpost[1,:] = tempup
			self.u_dispost = np.zeros(self.u_parpost.shape)
			self.u_dispost[1,:] = 1
		
	def do_sensitivity(self):
		if self.sensitype == 'k-s':
			self.s_result = []
			for ii in range(self.pvalue_a.shape[1]):
				d,p = ks_2samp(self.pvalue_a[:,ii],self.pvalue_r[:,ii])
				self.s_result.append([d,p])
				
		elif self.sensitype == 'mlr':
			self.s_result = []
			inls0 = self.__generate_parlist()
			outls0 = []
			for sim in self.simlist:
				x = sim.search_result(self.postin[0],self.postin[1],self.postin[2])				
				outls0.append(x)
			inls = np.insert(np.array(inls0),0,values=np.ones(self.simtime), axis=1)
			outls = np.array(outls0)
			ata = np.dot(inls.T,inls)
			atb = np.dot(inls.T,outls)
			coef = np.linalg.solve(ata,atb)
			sigma = (np.sum((np.dot(inls,coef) -outls) ** 2)/(inls.shape[0]-inls.shape[1])) ** 0.5
			cc = np.diag(np.linalg.inv(ata)) ** 0.5 * sigma
			tstat = coef/cc
			pvalue = (1 - t.cdf(tstat,(inls.shape[0]-inls.shape[1]))) * 2
			pvalue[pvalue>1] = 2 - pvalue[pvalue>1]
			for tup in zip(tstat,pvalue):
				self.s_result.append(tup)
			self.s_result.pop(0)
			
		elif self.sensitype == 'sobol':
			basen = int(self.simtime/(len(self.parnamels)+2))
			start = 0
			self.s_result = []
			outls0 = []
			for sim in self.simlist:
				x = sim.search_result(self.postin[0],self.postin[1],self.postin[2])				
				outls0.append(x)
			outls0 = np.array(outls0)
			vary = np.std(outls0[0:basen])
			for k in range(0,len(self.parnamels)):
				varx0 = np.sum(outls0[basen:2*basen]*(outls0[basen*(k+2):basen*(k+3)]-outls0[0:basen]))/basen 
				ex0 = np.sum((outls0[0:basen]-outls0[basen*(k+2):basen*(k+3)])**2)/(2*basen)			
				self.s_result.append([varx0/vary, ex0/vary])	
	def do_display(self):
		self.s_result = np.round(self.s_result,3)
		self.__draw_sensi_pic()
		self.__draw_uncer_parpic()
		if self.simlist[0].stalist[0].sim:
			self.__draw_uncer_varpic()
			
	def write_output(self):
		with open('./OUTPUT/CALI_OUT/TEXT/wholeoutput.csv','w') as ff:
			writer = csv.writer(ff, lineterminator='\n') 
			writer.writerow(self.__insert(list(self.parnamels), 'ParName'))
			writer.writerow(self.__insert([x[1] for x in self.parupdown], 'ParUp'))
			writer.writerow(self.__insert([x[0] for x in self.parupdown], 'ParDown'))
			writer.writerow(self.__insert(list(self.stanamels), 'StationName'))
			writer.writerow(self.__insert(list(self.varnamels), 'VarName'))
			writer.writerow(self.__insert(list(self.likedic), 'Function'))
			writer.writerow(self.__insert([a[0] for a in list(self.likedic.values())], 'Threshold'))
			writer.writerow(self.__insert([a[1] for a in list(self.likedic.values())], 'AcceptWay'))
			writer.writerow(['Simulations'])
			ls = self.__generate_title()
			a = self.__insert(self.parnamels, 'Number')
			a.extend(ls)
			a.append('Accept?')
			writer.writerow(a)
			ii = 1
			for sim in self.simlist:
				writer.writerow(self.__insert(sim.generate_result(),ii))
				ii += 1		
				
	def accept_reject(self):
		accls = [sim for sim in self.simlist if sim.accept == 'yes']
		rejls = [sim for sim in self.simlist if sim.accept == 'no']
		simnum_a = [sim.num for sim in accls]
		pvalue_a = self.__generate_parlist(accls)
		pvalue_r = self.__generate_parlist(rejls)
		
		like_a = []	
		base = -11060424
		sim0 = self.simlist[0]
		for sim in self.simlist:
			temp = sim.search_result(self.postin[0],self.postin[1],self.postin[2])
			if base == -11060424:
				base = temp
			if self.likedic[self.postin[2]][1] == '+':
				if temp > base:
					sim0 = sim
					base = temp
			elif self.likedic[self.postin[2]][1] == '-':
				if temp < base:
					sim0 = sim
					base = temp		
			elif self.likedic[self.postin[2]][1] == '[]':
				if np.abs(temp) < np.abs(base):
					sim0 = sim
					base = temp	
			if self.uncertype == 'sufi2':
				like_a.append(temp)
			elif self.uncertype == 'glue':
				if sim in accls:
					like_a.append(temp)
		self.best_sim = sim0
		best_par = [sim0.pardict[x] for x in self.parnamels]
		return simnum_a, pvalue_a, pvalue_r, like_a, best_par

	def save_pic(self,pic,name):
		try:
			pic.savefig('./OUTPUT/CALI_OUT/PIC/'+name, dpi=300, bbox_inches='tight')
			return 'done'
		except:
			return 'not'

	def save_sensi(self):
		try:
			with open('./OUTPUT/CALI_OUT/TEXT/s_result.csv','w') as ff:
				writer = csv.writer(ff, lineterminator='\n') 
				writer.writerow(['parname','stats value','p-value'])
				for ii in range(len(self.parnamels)):
					ls = [self.parnamels[ii], str(self.s_result[ii][0]), str(self.s_result[ii][1])]
					writer.writerow(ls)
			return 'done'
		except:
			return 'not'
			
	def save_parpost(self):
		try:	
			with open('./OUTPUT/CALI_OUT/TEXT/posterior.csv','w') as ff:
				writer = csv.writer(ff, lineterminator='\n') 
				for ii in range(len(self.parnamels)):
					if self.uncertype != 'sufi2':
						writer.writerow([str(self.parnamels[ii]),str(self.parupdown[ii][0]),'0'])
						for jj in range(self.u_parpost.shape[0]):
							writer.writerow(['',str(self.u_parpost[jj,ii]),str(self.u_dispost[jj,ii])])
						writer.writerow(['',str(self.parupdown[ii][1]),'1'])
					else:
						writer.writerow([str(self.parnamels[ii]),str(self.u_parpost[0,ii]),str(self.u_dispost[0,ii])])
						writer.writerow(['',str(self.u_parpost[1,ii]),str(self.u_dispost[1,ii])])
			with open('./OUTPUT/CALI_OUT/TEXT/acceptpar.txt','w') as ff:
				string = '\t'.join(self.parnamels)
				ff.write('NUM\t'+string+'\n')
				k = 1
				for line in self.pvalue_a:
					ff.write(str(k))
					for ii in line:
						ff.write('\t'+str(ii))
					ff.write('\n')
					k += 1
			with open('./OUTPUT/CALI_OUT/TEXT/bestpar.txt','w') as ff:
				string = '\t'.join(self.parnamels)
				ff.write('NUM\t'+string+'\n')
				ff.write('BEST')
				for ii in self.bestpar:
					ff.write('\t'+str(ii))
				ff.write('\n')		
			return 'done'
		except:
			return 'not'	
			
	def save_varpost(self):
		try:	
			with open('./OUTPUT/CALI_OUT/TEXT/95PPU.csv','w') as ff:
				writer = csv.writer(ff, lineterminator='\n') 
				writer.writerow(self.__insert(self.postin[0],'StaName'))
				writer.writerow(self.__insert(self.p_factor,'P-factor'))
				writer.writerow(self.__insert(self.r_factor,'R-factor'))			
			with open('./OUTPUT/CALI_OUT/TEXT/bestpar_calc.csv','w') as ff:
				writer = csv.writer(ff, lineterminator='\n')
				out = []
				ways = ['me','pbias','mae','mse','rmse','rsr','ns','r2','r-s','kge']
				for way in ways:
					out.append(self.__insert(self.stanamels,way))
					for varname in self.varnamels:
						temp = [varname]
						for sta in self.best_sim.stalist:
							try:
								sta.calc_do(way,varname,np.array(sta.sim[varname]),np.array(sta.obs[varname]))
								temp.append(sta.calc[varname][way])
							except:
								temp.append('')
						out.append(temp)
					out.append('')
				for line in out:
					writer.writerow(line)
			return 'done'
		except:
			return 'not'		
			
			
			
	def __generate_parlist(self,ls = []):
		back =[]
		if not ls:
			ls = self.simlist
		for sim in ls:
			back.append([sim.pardict[x] for x in self.parnamels])
		return back
	
	
	def __generate_title(self):
		a = self.simlist[0]
		ls = []
		for sta in a.stalist:
			for var in list(sta.calc):
				for func in list(sta.calc[var]):
					string = sta.staname + '+' + var + '+' + func
					ls.append(string)
		return ls
		
	def __insert(self, ls, value, place = 0):
		a = copy.deepcopy(ls)
		a.insert(place, value)
		return a
		
	def __draw_sensi_pic(self):	
		srarr = np.array(self.s_result)
		if self.sensitype == 'k-s':
			slb = 'k-s stats'
			vlb = 'p value'
		elif self.sensitype == 'mlr':
			slb = 't stats'
			vlb = 'p value'
		elif self.sensitype == 'sobol':
			slb = 'first order'
			vlb = 'global'
		if 8/len(self.parnamels) < 0.8:
			total_height = 8/len(self.parnamels)
		else: 
			total_height = 0.8
		n_bar = 2
		height = total_height/n_bar
		x = range(len(self.parnamels))
		

		plt.rcParams['figure.figsize'] = (5.0, 5.0)	
		plt.rcParams['figure.dpi'] = 100
		plt.rcParams['savefig.dpi'] = 300
		
		fig,ax1 = plt.subplots(1,1)
		ax1.set_xlim(-np.max(np.abs(srarr[:,0]))-1,np.max(np.abs(srarr[:,0]))+1)
		ax1.set_xlabel(slb)
		bar1 = ax1.barh(x, srarr[:,0], height=-height, fc='y')
		for i in range(len(x)):
			x[i] = x[i] + height
		ax2 = ax1.twiny()
		ax2.set_xlim(-1,1)
		ax2.set_xlabel(vlb)
		bar2 = ax2.barh(x, srarr[:,1], height=-height, fc='c',  tick_label=self.parnamels, align='center')

		ax1.legend([bar1,bar2],[slb,vlb], prop={'size': 12}) 
		ax1.spines['right'].set_visible(False)  
		ax2.spines['right'].set_visible(False) 		
		ax1.spines['left'].set_position(('data',0))
		ax2.spines['left'].set_position(('data',0)) 	
		self.sensi_pic = plt.gcf()
		plt.close(self.sensi_pic)	
		
	def __draw_uncer_parpic(self):
		inls = np.array(self.__generate_parlist())
		outls = []
		for sim in self.simlist:
			x = sim.search_result(self.postin[0],self.postin[1],self.postin[2])				
			outls.append(x)
		outls = np.array(outls)
		self.par_plot_picls = []
		for ii in range(len(self.parnamels)):
			a = plt.figure()		
			ax = plt.gca()
			plt.scatter(inls[:,ii],outls)
			plt.title(self.parnamels[ii])
			plt.legend(self.parnamels[ii])	
			
			xticks = ax.get_xticks()
			yticks = ax.get_yticks()
			xticks[0] = self.parupdown[ii][0]
			xticks[-1] = self.parupdown[ii][1]
			ax.set_xticks(xticks)
			ax.set_xlim(self.parupdown[ii][0],self.parupdown[ii][1])
			ax.set_ylim(yticks[0],yticks[-1])		

			self.par_plot_picls.append(a)
			plt.close(a)
		if self.pvalue_a.shape[0] > 50:
			histnum = int(inls.shape[0]/10)
		else:
			histnum = 5
		self.posterior_picls = []
		for ii in range(len(self.parnamels)):
			a = plt.figure()
			plt.hist(self.pvalue_a[:,ii],histnum,range=(self.parupdown[ii][0],self.parupdown[ii][1]),edgecolor='black',facecolor='green',alpha=0.5)
			plt.title(self.parnamels[ii])
			
			ax = plt.gca()
			xticks = [self.parupdown[ii][0] + kk/histnum * (self.parupdown[ii][1]-self.parupdown[ii][0]) for kk in range(histnum+1) ]
			ax.set_xticks(xticks)
			ax.set_xlim(self.parupdown[ii][0],self.parupdown[ii][1])	

			self.posterior_picls.append(a)
			plt.close(a)
	
	def __draw_uncer_varpic(self):
		if self.simlist[0].stalist[0].obs:
			plt.rcParams['figure.figsize'] = (5.0, 2.0)	
			self.time_pic = []
			self.time_best = []
			self.p_factor = []
			self.r_factor = []
			simlist = [s for s in self.simlist if s.num in self.simnum_a]
			for name in self.postin[0]:
				obs = []
				sim = []
				bestsim = []
				for s in simlist:
					for sta in s.stalist:
						if sta.staname == name:
							sim.append(sta.sim[self.postin[1]])
							if obs == []:
								obs = sta.obs[self.postin[1]]
							if s == self.best_sim:
								bestsim = sta.sim[self.postin[1]]
				sim = np.array(sim)
				max = np.percentile(sim,97.5,axis=0)
				min = np.percentile(sim,2.5,axis=0)
				xaxis = range(1,len(max)+1)
				p = [x for x in np.all(np.array([obs>=min,obs<=max]),axis=0) if x]
				self.p_factor.append(float(np.round(len(p)/len(obs),3)))
				r = np.mean(max-min)/np.std(obs)
				self.r_factor.append(float(np.round(r,3)))			
				a = plt.figure()
				plt.fill_between(xaxis,min,max,color='cyan', alpha=0.75)
				plt.plot(xaxis,obs,'r.-')
				plt.title(name)
				ax = plt.gca()
				ax.legend(['obs','95PPU'], prop={'size': 12}) 
				ax.spines['top'].set_visible(False) 
				ax.spines['right'].set_visible(False)  
				self.time_pic.append(a)
				plt.close(a)			
				b = plt.figure()
				plt.plot(xaxis,obs,'r-')
				plt.plot(xaxis,bestsim,'g-')
				ax = plt.gca()
				plt.title(name)
				ax.legend(['obs','bestsim'], prop={'size': 12}) 
				ax.spines['top'].set_visible(False) 
				ax.spines['right'].set_visible(False)  	
				self.time_best.append(b)
				plt.close(b)
