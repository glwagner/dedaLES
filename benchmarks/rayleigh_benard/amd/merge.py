import numpy as np
import matplotlib.pyplot as plt
import h5py
from dedalus import public as de
from dedalus.extras import flow_tools
import time
from IPython import display
import sys

from dedalus.tools import post
print('merging files from')
print(sys.argv[1])
post.merge_process_files(sys.argv[1], cleanup=True)
#when merging the file name and subfiles must share a common name
#in other words one cannot rename the top level directory into something else and expect it to work

