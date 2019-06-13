from langmuir import *
from frmt import *

t = get_table('darian-marholm sphere')

top = t['axes'][2] # R
left = t['axes'][3] # V
vals = t['values'][0][0].T

t = np.zeros((vals.shape[0]+1, vals.shape[1]+1))
t[1:,1:] = vals
t[0, 1:] = top
t[1:, 0] = left

print_table(t, '>', [['V\R','{:.1f}'],['{:.0f}','{:.3f}']])
