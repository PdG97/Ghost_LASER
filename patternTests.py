# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 11:40:22 2021

@author: Labor_NLO-Admin
"""
import numpy as np
import Pattern as pt
p = pt.PatternHelper()
pattern = p.getFreqSweepLayout(430, 700, 10)
pattern.plot()
