import sys
import numpy as np
import builtins
import os
import sys
from types import ModuleType, FunctionType
from gc import get_referents
from pympler.asizeof import asizeof as msize
from time import time as t


def write_variables():  
    cwd = os.getcwd()
    os.chdir('/home/mikkola/Documents/DeprojectionProject')
    variables = list(builtins.lv.items())
    l_vars = []
    for var, obj in variables:
        l_vars += [(var.ljust(30), obj/1e9)]
    l_vars = sorted(l_vars, key=lambda mem: mem[1], reverse=True)  
    with open('variable_memory.log', 'w+') as vlog:
        vlog.write('Variable name'.ljust(30) + 'Size'.ljust(15) + 'GB\n')
        for line in l_vars:
            vlog.write(line[0] + str(line[1]).ljust(15) + 'GB\n')

    os.chdir(cwd)
    return


def mod_dict(vdict):
    for var in vdict.keys():
        vdict[var] = getsize(vdict[var])
    return vdict

a = np.random.normal(size=(200,100,100))
#func()


# Custom objects know their class.
# Function objects seem to know way too much, including modules.
# Exclude modules as well.
BLACKLIST = type, ModuleType, FunctionType


def getsize(obj):
    """sum size of object & members."""
    if isinstance(obj, BLACKLIST):
        return 0
    seen_ids = set()
    size = 0
    objects = [obj]
    while objects:
        need_referents = []
        for obj in objects:
            if not isinstance(obj, BLACKLIST) and id(obj) not in seen_ids:
                seen_ids.add(id(obj))
                size += sys.getsizeof(obj)
                need_referents.append(obj)
        objects = get_referents(*need_referents)
    return size


builtins.lv = dict(locals())
print('locals(): ', getsize(builtins.lv)/1e9)
test = np.random.normal(size=(2000,100,1000))
print('test = ...: ', getsize(test)/1e9)
builtins.lv.update(mod_dict(dict(locals())))
print('locals(): ', getsize(builtins.lv)/1e9)
print(builtins.lv['test'])
write_variables()

