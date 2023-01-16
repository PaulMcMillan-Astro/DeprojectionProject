#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 21:05:44 2020

@author: mikkola
"""
import sys
import numpy as np
import curses
from decimal import Decimal
from curses import wrapper
import time

def map(stdscr, st,y, it_type):
        yvals = y
        ysort = np.argsort(y)
        y = np.concatenate([np.where(np.sort(y) == y[i])[0] for i in range(4)])
        if it_type == 'gss':
            ys = [""," "*6," "*12," "*18]
            for i in range(4):
                if i == 0:
                    ys[i] += "O-----|"
                elif i == 3:
                    ys[i] += "|-----O"
                else:
                    ys[i] += "|--O--|"
                    
            for k in range(3):
                ln = [6,12,18]
                list_i = [y[k],y[k+1]]
                for i in np.arange(np.min(list_i)+1,np.max(list_i),1):
                    temp_s = list(ys[ysort[i]].ljust(30))
                    temp_s[ln[k]] = "|"
                    ys[ysort[i]] = "".join(temp_s) 
        elif it_type == 'ternary':
            ys = [""," "*4," "*12," "*20]
            for i in range(4):
                if i == 0:
                    ys[i] += "O---|"
                elif i == 3:
                    ys[i] += "|---O"
                else:
                    ys[i] += "|---O---|"
                    
            for k in range(3):
                ln = [4,12,20]
                list_i = [y[k],y[k+1]]
                for i in np.arange(np.min(list_i)+1,np.max(list_i),1):
                    temp_s = list(ys[ysort[i]].ljust(30))
                    temp_s[ln[k]] = "|"
                    ys[ysort[i]] = "".join(temp_s)             
        ln = 0
        for i in reversed(ysort):
            stdscr.addstr(st+ln, 8, ys[i], curses.color_pair(0) | curses.A_BOLD)
            ln += 1
        cn = 0    
        for i in range(4):  
            stdscr.addstr(st+ln+2, 0+cn, str("%.4e" % Decimal(yvals[i])).ljust(6), curses.color_pair(3))
            cn += 11
            
        

def gssstrings(stdscr):
    
    st = 51
    stdscr.addstr(0, 1, "* Using opt_alpha_gss\n", curses.color_pair(7) | curses.A_BOLD | curses.A_UNDERLINE)
    stdscr.addstr(2, st,("Types of steps:"), curses.color_pair(0))
    stdscr.addstr(3, st,("a = current range value, b = new range value"), curses.color_pair(0))
    #RIGHT
    st = 69
    stdscr.addstr(5,st-18,"Right", curses.color_pair(0) | curses.A_BOLD )
    stdscr.addstr(6,st,"a1 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(6,st+3,"--------- ", curses.color_pair(0))
    stdscr.addstr(6,st+13,"a2 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(6,st+16,"-- ", curses.color_pair(0))
    stdscr.addstr(6,st+19,"a3 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(6,st+22,"--------- ", curses.color_pair(0))
    stdscr.addstr(6,st+32,"a4", curses.color_pair(7) | curses.A_BOLD ) 
    stdscr.addstr(7,st+13,"||" +  "\\\\".rjust(6) +  "||".rjust(13), curses.color_pair(0))
    st = 79
    stdscr.addstr(8,st,"-> ", curses.color_pair(0))
    stdscr.addstr(8,st+3,"a2 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(8,st+6,"--- ", curses.color_pair(0))
    stdscr.addstr(8,st+10,"a3 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(8,st+13,"- ", curses.color_pair(0))
    stdscr.addstr(8,st+15,"b3 ", curses.color_pair(2) | curses.A_BOLD )
    stdscr.addstr(8,st+18,"--- ", curses.color_pair(0))
    stdscr.addstr(8,st+22,"a4", curses.color_pair(7) | curses.A_BOLD )
    #LEFT
    st = 69
    stdscr.addstr(10,st-18,"Left", curses.color_pair(0) | curses.A_BOLD )
    stdscr.addstr(11,st,"a1 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(11,st+3,"--------- ", curses.color_pair(0))
    stdscr.addstr(11,st+13,"a2 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(11,st+16,"-- ", curses.color_pair(0))
    stdscr.addstr(11,st+19,"a3 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(11,st+22,"--------- ", curses.color_pair(0))
    stdscr.addstr(11,st+32,"a4", curses.color_pair(7) | curses.A_BOLD )  
    stdscr.addstr(12,st,"||" +  "//".rjust(13)+  "||".rjust(6), curses.color_pair(0))
    st = 66
    stdscr.addstr(13,st,"-> ", curses.color_pair(0))
    stdscr.addstr(13,st+3,"a1 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(13,st+6,"--- ", curses.color_pair(0))
    stdscr.addstr(13,st+10,"b2 ", curses.color_pair(2) | curses.A_BOLD )
    stdscr.addstr(13,st+13,"- ", curses.color_pair(0))
    stdscr.addstr(13,st+15,"a2 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(13,st+18,"--- ", curses.color_pair(0))
    stdscr.addstr(13,st+22,"a3", curses.color_pair(7) | curses.A_BOLD )
    #UP
    st = 69
    stdscr.addstr(15,st-18,"Up", curses.color_pair(0) | curses.A_BOLD )    
    stdscr.addstr(16,st,"a1 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(16,st+3,"--------- ", curses.color_pair(0))
    stdscr.addstr(16,st+13,"a2 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(16,st+16,"-- ", curses.color_pair(0))
    stdscr.addstr(16,st+19,"a3 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(16,st+22,"--------- ", curses.color_pair(0))
    stdscr.addstr(16,st+32,"a4", curses.color_pair(7) | curses.A_BOLD ) 
    stdscr.addstr(16,st+35,"     ...      v", curses.color_pair(0))
    stdscr.addstr(17,st+19,("||" + "v".rjust(29)), curses.color_pair(0))
    st = 85
    stdscr.addstr(18,st,"-> ", curses.color_pair(0))
    stdscr.addstr(18,st+3,"a3 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(18,st+6,"------ ", curses.color_pair(0))
    stdscr.addstr(18,st+14,"b2 ", curses.color_pair(2) | curses.A_BOLD )
    stdscr.addstr(18,st+17,"--- ", curses.color_pair(0))
    stdscr.addstr(18,st+21,"b3 ", curses.color_pair(2) | curses.A_BOLD )
    stdscr.addstr(18,st+24,"------ (", curses.color_pair(0))
    stdscr.addstr(18,st+32,"a4" , curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(18,st+34,"+1)", curses.color_pair(0))
    #DOWN
    st = 52 
    stdscr.addstr(20,st-1,"Down", curses.color_pair(0) | curses.A_BOLD )
    stdscr.addstr(21,st,"v      ...       ", curses.color_pair(0))
    stdscr.addstr(21,st+17,"a1 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(21,st+20,"--------- ", curses.color_pair(0))
    stdscr.addstr(21,st+30,"a2 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(21,st+33,"-- ", curses.color_pair(0))
    stdscr.addstr(21,st+36,"a3 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(21,st+39,"--------- ", curses.color_pair(0))
    stdscr.addstr(21,st+49,"a4", curses.color_pair(7) | curses.A_BOLD ) 
    stdscr.addstr(22,st,("v" + "||".rjust(31)), curses.color_pair(0))
    st = 47
    stdscr.addstr(23,st,"->(", curses.color_pair(0))
    stdscr.addstr(23,st+4,"a1", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(23,st+6,"-1) ", curses.color_pair(0))
    stdscr.addstr(23,st+10,"------ ", curses.color_pair(0))
    stdscr.addstr(23,st+17,"b2 ", curses.color_pair(2) | curses.A_BOLD )
    stdscr.addstr(23,st+20,"--- ", curses.color_pair(0))
    stdscr.addstr(23,st+24,"b3 ", curses.color_pair(2) | curses.A_BOLD )
    stdscr.addstr(23,st+28,"------ ", curses.color_pair(0))
    stdscr.addstr(23,st+35,"a2", curses.color_pair(7) | curses.A_BOLD )

    stdscr.addstr(11, 7, ("a1".ljust(28) + "a4"), curses.color_pair(3) | curses.A_BOLD)
    stdscr.addstr(13, 7, ("b1".ljust(28) + "b4"), curses.color_pair(3) | curses.A_BOLD)
    return

def ternstrings(stdscr):
    
    st = 51
    stdscr.addstr(0, 1, "* Using opt_alpha_ternary", curses.color_pair(7) | curses.A_BOLD | curses.A_UNDERLINE)
    stdscr.addstr(2,st,("Types of steps:"), curses.color_pair(0))
    stdscr.addstr(3,st,("a = current range value, b = new range value"), curses.color_pair(0))
    #RIGHT
    st = 69
    stdscr.addstr(5,st-18,"Right", curses.color_pair(0) | curses.A_BOLD )
    stdscr.addstr(6,st,"a1 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(6,st+3,"----- ", curses.color_pair(0))
    stdscr.addstr(6,st+9,"a2 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(6,st+12,"----- ", curses.color_pair(0))
    stdscr.addstr(6,st+18,"a3 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(6,st+21,"----- ", curses.color_pair(0))
    stdscr.addstr(6,st+27,"a4", curses.color_pair(7) | curses.A_BOLD ) 
    stdscr.addstr(7,st+9,("||" +  "||".rjust(18)), curses.color_pair(0))
    st = 75
    stdscr.addstr(8,st,"-> ", curses.color_pair(0))
    stdscr.addstr(8,st+3,"a2 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(8,st+6,"-- ", curses.color_pair(0))
    stdscr.addstr(8,st+9,"b2 ", curses.color_pair(2) | curses.A_BOLD )
    stdscr.addstr(8,st+12,"-- ", curses.color_pair(0))
    stdscr.addstr(8,st+15,"b3 ", curses.color_pair(2) | curses.A_BOLD )
    stdscr.addstr(8,st+18,"-- ", curses.color_pair(0))
    stdscr.addstr(8,st+21,"a4", curses.color_pair(7) | curses.A_BOLD )
    #LEFT
    st = 69
    stdscr.addstr(10,st-18,"Left", curses.color_pair(0) | curses.A_BOLD )
    stdscr.addstr(11,st,"a1 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(11,st+3,"----- ", curses.color_pair(0))
    stdscr.addstr(11,st+9,"a2 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(11,st+12,"----- ", curses.color_pair(0))
    stdscr.addstr(11,st+18,"a3 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(11,st+21,"----- ", curses.color_pair(0))
    stdscr.addstr(11,st+27,"a4", curses.color_pair(7) | curses.A_BOLD ) 
    stdscr.addstr(12,st,("||" + "||".rjust(18)), curses.color_pair(0))
    st = 66
    stdscr.addstr(13,st,"-> ", curses.color_pair(0))
    stdscr.addstr(13,st+3,"a1 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(13,st+6,"-- ", curses.color_pair(0))
    stdscr.addstr(13,st+9,"b2 ", curses.color_pair(2) | curses.A_BOLD )
    stdscr.addstr(13,st+12,"-- ", curses.color_pair(0))
    stdscr.addstr(13,st+15,"b3 ", curses.color_pair(2) | curses.A_BOLD )
    stdscr.addstr(13,st+18,"-- ", curses.color_pair(0))
    stdscr.addstr(13,st+21,"a3", curses.color_pair(7) | curses.A_BOLD )
    #UP
    st = 69
    stdscr.addstr(15,st-18,"Up", curses.color_pair(0) | curses.A_BOLD )    
    stdscr.addstr(16,st,"a1 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(16,st+3,"----- ", curses.color_pair(0))
    stdscr.addstr(16,st+9,"a2 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(16,st+12,"----- ", curses.color_pair(0))
    stdscr.addstr(16,st+18,"a3 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(16,st+21,"----- ", curses.color_pair(0))
    stdscr.addstr(16,st+27,"a4", curses.color_pair(7) | curses.A_BOLD ) 
    stdscr.addstr(16,st+29,"   ...   v", curses.color_pair(0))
    stdscr.addstr(17,st+18,("||" + "v".rjust(19)), curses.color_pair(0))
    st = 84
    stdscr.addstr(18,st,"-> ", curses.color_pair(0))
    stdscr.addstr(18,st+3,"a3 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(18,st+6,"-- ", curses.color_pair(0))
    stdscr.addstr(18,st+9,"b2 ", curses.color_pair(2) | curses.A_BOLD )
    stdscr.addstr(18,st+12,"-- ", curses.color_pair(0))
    stdscr.addstr(18,st+15,"b3 ", curses.color_pair(2) | curses.A_BOLD )
    stdscr.addstr(18,st+18,"-- (", curses.color_pair(0))
    stdscr.addstr(18,st+22,"a4" , curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(18,st+24,"+1)", curses.color_pair(0))
    #DOWN
    st = 52
    stdscr.addstr(20,st-1,"Down", curses.color_pair(0) | curses.A_BOLD )
    stdscr.addstr(21,st,"v   ...   ", curses.color_pair(0))
    stdscr.addstr(21,st+10,"a1 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(21,st+13,"----- ", curses.color_pair(0))
    stdscr.addstr(21,st+19,"a2 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(21,st+22,"----- ", curses.color_pair(0))
    stdscr.addstr(21,st+28,"a3 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(21,st+31,"----- ", curses.color_pair(0))
    stdscr.addstr(21,st+37,"a4", curses.color_pair(7) | curses.A_BOLD )  
    stdscr.addstr(22,st,("v" + "||".rjust(20)), curses.color_pair(0))
    st = 46
    stdscr.addstr(23,st,"-> (", curses.color_pair(0))
    stdscr.addstr(23,st+4,"a1", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(23,st+6,"-1) ", curses.color_pair(0))
    stdscr.addstr(23,st+10,"-- ", curses.color_pair(0))
    stdscr.addstr(23,st+13,"b2 ", curses.color_pair(2) | curses.A_BOLD )
    stdscr.addstr(23,st+16,"-- ", curses.color_pair(0))
    stdscr.addstr(23,st+19,"b3 ", curses.color_pair(2) | curses.A_BOLD )
    stdscr.addstr(23,st+22,"-- ", curses.color_pair(0))
    stdscr.addstr(23,st+25,"a2", curses.color_pair(7) | curses.A_BOLD )

    stdscr.addstr(11, 7, ("a1".ljust(28) + "a4"), curses.color_pair(3) | curses.A_BOLD)
    stdscr.addstr(13, 7, ("b1".ljust(28) + "b4"), curses.color_pair(3) | curses.A_BOLD)
    return

def tenstepstrings(stdscr):
    st = 51
    stdscr.addstr(0, 1, "* Using opt_alpha", curses.color_pair(7) | curses.A_BOLD | curses.A_UNDERLINE)
    stdscr.addstr(2,st,("Types of steps:"), curses.color_pair(0))
    stdscr.addstr(3,st,("a = current range value, b = new range value"), curses.color_pair(0))
    stdscr.addstr(4,st,("Magenta = a with lowest MISE"), curses.color_pair(0))
    #ZOOM
    st = 64
    stdscr.addstr(5,st-13,"Zoom", curses.color_pair(0) | curses.A_BOLD )
    for n in range(1,10):
        stdscr.addstr(6,st,str("a%s " % n), curses.color_pair(7) | curses.A_BOLD ); st += 3
        stdscr.addstr(6,st,"- ", curses.color_pair(0)); st += 2
    stdscr.addstr(6,st,str("a%s " % 10), curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(6,st-25,str("a%s " % 5), curses.color_pair(6))
    stdscr.addstr(7,st-28,("/" +  "\\".rjust(7)), curses.color_pair(0))
    st = 77
    stdscr.addstr(8,st,"b1 ", curses.color_pair(2) | curses.A_BOLD )
    stdscr.addstr(8,st+3,"<-..10x..-> ", curses.color_pair(0))
    stdscr.addstr(8,st+15,"b10", curses.color_pair(2) | curses.A_BOLD )
    
    #UP
    st = 64
    stdscr.addstr(10,st-13,"Up", curses.color_pair(0) | curses.A_BOLD )
    for n in range(1,10):
        stdscr.addstr(11,st,str("a%s " % n), curses.color_pair(7) | curses.A_BOLD ); st += 3
        stdscr.addstr(11,st,"- ", curses.color_pair(0)); st += 2
    stdscr.addstr(11,st,str("a%s " % 10), curses.color_pair(6))
    stdscr.addstr(11,st+3,"-------", curses.color_pair(0))
    stdscr.addstr(12,st-5,("||" +  "\\".rjust(14)), curses.color_pair(0))
    st = 104
    stdscr.addstr(13,st,"a9 ", curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(13,st+3,"<-...10x...-> ", curses.color_pair(0))
    stdscr.addstr(13,st+17,"b10", curses.color_pair(2) | curses.A_BOLD )
    #DOWN
    st = 64
    stdscr.addstr(15,st-13,"Down", curses.color_pair(0) | curses.A_BOLD )
    for n in range(1,10):
        stdscr.addstr(16,st,str("a%s " % n), curses.color_pair(7) | curses.A_BOLD ); st += 3
        stdscr.addstr(16,st,"- ", curses.color_pair(0)); st += 2
    stdscr.addstr(16,st,str("a%s " % 10), curses.color_pair(7) | curses.A_BOLD )
    stdscr.addstr(16,st-45,str("a%s " % 1), curses.color_pair(6))
    stdscr.addstr(16,st-53,"--------", curses.color_pair(0))
    stdscr.addstr(17,st-54,("/" +  "||".rjust(15)), curses.color_pair(0))
    st = 52
    stdscr.addstr(18,st,"b1 ", curses.color_pair(2) | curses.A_BOLD )
    stdscr.addstr(18,st+3,"<-...10x...-> ", curses.color_pair(0))
    stdscr.addstr(18,st+17,"a2", curses.color_pair(7) | curses.A_BOLD )

    stdscr.addstr(11, 6, ("a1".ljust(29) + "a10"), curses.color_pair(3) | curses.A_BOLD)
    stdscr.addstr(13, 6, ("b1".ljust(29) + "b10"), curses.color_pair(3) | curses.A_BOLD)
    return


def make_string(stdscr, logalpha0, opt_it, mise_it, mise_it_tot, mise_diff, xrange, xrange_old, yrange_old, t, step, it_type):
    x1,x4,x1new,x4new = xrange_old[0],xrange_old[-1],xrange[0],xrange[-1]
    
    stdscr.erase(); stdscr.refresh()
    stdscr.addstr(1, 1, str('logalpha0    : %.2f' % logalpha0), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(2, 1, str('opt it       : %s' % opt_it), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(3, 1, str('mise it      : %s (%s tot)' % (mise_it,mise_it_tot)), curses.color_pair(0) | curses.A_BOLD) 
    stdscr.addstr(4, 1, str('mise diff    : %s' % (np.around(abs(x4new-x1new),2))), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(5, 1, str('time elapsed : %s hrs' % (np.around(t,3))), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(10, 1, (str("%.6f"%x1).ljust(10).ljust(28) + str("%.6f"%x4)), curses.color_pair(3))
    stdscr.addstr(12, 13, ("\/".ljust(7) + step.ljust(5) + "\/".rjust(5)), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(14, 1, (str("%.6f"%x1new).ljust(10).ljust(28) + str("%.6f"%x4new)), curses.color_pair(3))
    
    if it_type == 'gss':
        gssstrings(stdscr)
        st = 20
        stdscr.addstr(st-2, 1, 'MISE structure:', curses.color_pair(0) | curses.A_BOLD)
        map(stdscr, st,yrange_old, it_type)
    elif it_type == 'ternary':
        ternstrings(stdscr)
        st = 20
        stdscr.addstr(st-2, 1, 'MISE structure:', curses.color_pair(0) | curses.A_BOLD)
        map(stdscr, st,yrange_old, it_type)
    elif it_type == 'tenstep':
        tenstepstrings(stdscr)
        
    
    for i in range(25):
        stdscr.addstr(i,44,"|", curses.color_pair(0))
    stdscr.refresh()    

def debugger():
    curses.start_color()
    curses.use_default_colors()
    for i in range(0, curses.COLORS):
        curses.init_pair(i + 1, i, -1)



