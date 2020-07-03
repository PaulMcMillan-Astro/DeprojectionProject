from sys import argv
import os
os.chdir('/home/mikkola/Documents/DeprojectionProject')
from Deproject_v1_0 import *
import time
from astropy.io.ascii import read as tableread
from datetime import date
import builtins
import string
from termcolor import colored
from Deproject_test import sanity_check
from Deproject_plots import plot_fv,plot_L

#import sys
#import numpy as np
#import curses
#from curses import wrapper
#import time

#def map(st,y):
   #ysort = np.argsort(y)
   #y = np.concatenate([np.where(np.sort(y) == y[i])[0] for i in range(4)])
   #ys = [""," "*4," "*12," "*20]
   #for i in range(4):
       #if i == 0:
           #ys[i] += "O---|"
       #elif i == 3:
           #ys[i] += "|---O"
       #else:
           #ys[i] += "|---O---|"
           
   #for k in range(3):
       #ln = [4,12,20]
       #list_i = [y[k],y[k+1]]
       #for i in np.arange(np.min(list_i)+1,np.max(list_i),1):
           #temp_s = list(ys[ysort[i]].ljust(30))
           #temp_s[ln[k]] = "|"
           #ys[ysort[i]] = "".join(temp_s) 
   #ln = 0
   #for i in reversed(ysort):
       #stdscr.addstr(st+ln, 8, ys[i], curses.color_pair(0) | curses.A_BOLD)
       #ln += 1

#def f(stdscr):
   #logalpha0, opt_it, mise_it, mise_it_tot, mise_diff, xrange, xrange_old, yrange_old, t, step = 1.222, 1, 1, 1, 1, [10,7,8,15], [1,2,3,4], [10,19,8,15], 200, 'Down' 
   #x1,x4,x1new,x4new = xrange_old[0],xrange_old[-1],xrange[0],xrange[-1]
   
   #st = 20
   #stdscr.addstr(st-2, 1, 'MISE structure:', curses.color_pair(0) | curses.A_BOLD)
   #map(st,yrange_old)
   
   #st = 51
   #stdscr.addstr(0, 1, "* Using opt_alpha_ternary", curses.color_pair(7) | curses.A_BOLD | curses.A_UNDERLINE)
   #stdscr.addstr(2,st,("Types of steps:"), curses.color_pair(0))
   #stdscr.addstr(3,st,("a = current range value, b = new range value"), curses.color_pair(0))
   #stdscr.addstr(4,st,("Magenta = a with lowest MISE"), curses.color_pair(0))
   #ZOOM
   #st = 64
   #stdscr.addstr(5,st-13,"Zoom", curses.color_pair(0) | curses.A_BOLD )
   #for n in range(1,10):
       #stdscr.addstr(6,st,str("a%s " % n), curses.color_pair(7) | curses.A_BOLD ); st += 3
       #stdscr.addstr(6,st,"- ", curses.color_pair(0)); st += 2
   #stdscr.addstr(6,st,str("a%s " % 10), curses.color_pair(7) | curses.A_BOLD )
   #stdscr.addstr(6,st-25,str("a%s " % 5), curses.color_pair(6))
   #stdscr.addstr(7,st-28,("/" +  "\\".rjust(7)), curses.color_pair(0))
   #st = 77
   #stdscr.addstr(8,st,"b1 ", curses.color_pair(2) | curses.A_BOLD )
   #stdscr.addstr(8,st+3,"<-..10x..-> ", curses.color_pair(0))
   #stdscr.addstr(8,st+15,"b10", curses.color_pair(2) | curses.A_BOLD )
   
   #UP
   #st = 64
   #stdscr.addstr(10,st-13,"Up", curses.color_pair(0) | curses.A_BOLD )
   #for n in range(1,10):
       #stdscr.addstr(11,st,str("a%s " % n), curses.color_pair(7) | curses.A_BOLD ); st += 3
       #stdscr.addstr(11,st,"- ", curses.color_pair(0)); st += 2
   #stdscr.addstr(11,st,str("a%s " % 10), curses.color_pair(6))
   #stdscr.addstr(11,st+3,"-------", curses.color_pair(0))
   #stdscr.addstr(12,st-5,("||" +  "\\".rjust(14)), curses.color_pair(0))
   #st = 104
   #stdscr.addstr(13,st,"a9 ", curses.color_pair(7) | curses.A_BOLD )
   #stdscr.addstr(13,st+3,"<-...10x...-> ", curses.color_pair(0))
   #stdscr.addstr(13,st+17,"b10", curses.color_pair(2) | curses.A_BOLD )
   #DOWN
   #st = 64
   #stdscr.addstr(15,st-13,"Down", curses.color_pair(0) | curses.A_BOLD )
   #for n in range(1,10):
       #stdscr.addstr(16,st,str("a%s " % n), curses.color_pair(7) | curses.A_BOLD ); st += 3
       #stdscr.addstr(16,st,"- ", curses.color_pair(0)); st += 2
   #stdscr.addstr(16,st,str("a%s " % 10), curses.color_pair(7) | curses.A_BOLD )
   #stdscr.addstr(16,st-45,str("a%s " % 1), curses.color_pair(6))
   #stdscr.addstr(16,st-53,"--------", curses.color_pair(0))
   #stdscr.addstr(17,st-54,("/" +  "||".rjust(15)), curses.color_pair(0))
   #st = 52
   #stdscr.addstr(18,st,"b1 ", curses.color_pair(2) | curses.A_BOLD )
   #stdscr.addstr(18,st+3,"<-...10x...-> ", curses.color_pair(0))
   #stdscr.addstr(18,st+17,"a2", curses.color_pair(7) | curses.A_BOLD )

   #for i in range(25):
       #stdscr.addstr(i,44,"|", curses.color_pair(0))

   #stdscr.addstr(i,44,"|", curses.color_pair(0))
   #stdscr.refresh() 

   #stdscr.getkey() 
   
   #stdscr.erase()
   #stdscr.addstr(6,1,"Press any key 3 times to log output and close program!", curses.color_pair(1))
   #stdscr.addstr(7,1,"Press any key 3 times to log output and close program!", curses.color_pair(2) | curses.A_BOLD)
   #stdscr.addstr(8,1,"Press any key 3 times to log output and close program!", curses.color_pair(3) | curses.A_BLINK)
   #stdscr.addstr(9,1,"Press any key 3 times to log output and close program!", curses.color_pair(4) | curses.A_BLINK)
   #stdscr.addstr(10,1,"Press any key 3 times to log output and close program!", curses.color_pair(5) | curses.A_REVERSE)
   #stdscr.addstr(11,1,"Press any key 3 times to log output and close program!", curses.color_pair(6) | curses.A_STANDOUT)
   #stdscr.addstr(12,1,"Press any key 3 times to log output and close program!", curses.color_pair(8) | curses.A_UNDERLINE)
   #stdscr.refresh() 
   #stdscr.getkey()   
   #stdscr.getkey() 
   #stdscr.getkey()
 
#stdscr = curses.initscr()
#curses.start_color()
#curses.use_default_colors()
#for i in range(0, curses.COLORS):
   #curses.init_pair(i + 1, i, -1)

#wrapper(f)

#import curses
#from curses import wrapper
#stdscr = curses.initscr()

#wrapper(f)
#try:
#   f()
#    time.sleep(1)
#except:
#    curses.echo()
#    curses.nocbreak()
#    curses.endwin()
#    raise TypeError("An error")
#else:
#    curses.echo()
#    curses.nocbreak()
#    curses.endwin()
