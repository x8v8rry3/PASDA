# -*- coding: utf-8 -*-
from __future__ import division	
import re

def std_split(string, filt = True, tofloat = False, cond = r'\s+'):
	ls = re.split(cond,string)
	if filt:
		ls = filter(lambda x:x,ls)
	if tofloat:
		ls = [[float(x),x][x == ''] for x in ls]	
	return ls
