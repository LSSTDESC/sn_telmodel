from .version import __version__
import os
import matplotlib.pyplot as plt

# requested infos for throughputs
"""
throughputs_dir = os.path.join(os.getenv("PWD"), "throughputs")
if not os.path.isdir(throughputs_dir):
    cmd = 'git clone https://github.com/lsst/throughputs'
    os.system(cmd)

os.environ['LSST_THROUGHPUTS_BASELINE'] = '{}/{}'.format(
    throughputs_dir, 'baseline')
os.environ['THROUGHPUTS_DIR'] = throughputs_dir

print('Reading throughputs from',
      os.environ['LSST_THROUGHPUTS_BASELINE'], os.environ['THROUGHPUTS_DIR'])
"""
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['figure.titlesize'] = 20
plt.rcParams['legend.fontsize'] = 20
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['figure.titleweight'] = 'bold'
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 20
