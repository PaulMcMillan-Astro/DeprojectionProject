from Deproject_v1_0 import KvalsBlockMethod, KvalsNumpyMethod
import pickle
import numpy as np

with open('vars.txt', 'rb') as f:
    var = pickle.loads(f.read()) # Keys: var['pvals'], var['rhat'], var['vmin'], var['dv'], var['n']
with open('output.txt', 'rb') as f:
    old_out = pickle.loads(f.read()) # Keys: output['numpy'], output['block'] for the two functions

# Convert from sparse matrices to arrays
old_out.update((label, item.toarray()) for label, item in old_out.items())

# Save new output, might as well make it a dict too..
new_out = {'numpy' : KvalsNumpyMethod(var['pvals'], var['rhat'], var['vmin'], var['dv'], var['n']).toarray(),
           'block' : KvalsBlockMethod(var['pvals'], var['rhat'], var['vmin'], var['dv'], var['n'], var['Nblock']).toarray()
            }


def test_kvals():
    # Numpy tests
    assert type(old_out['numpy']) == type(new_out['numpy'])
    assert old_out['numpy'].shape == new_out['numpy'].shape
    assert old_out['numpy'].sum() == new_out['numpy'].sum()
    
    # Block tests
    assert type(old_out['block']) == type(new_out['block'])
    assert old_out['block'].shape == new_out['block'].shape
    assert old_out['block'].sum() == new_out['block'].sum()