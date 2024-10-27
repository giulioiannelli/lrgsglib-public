#
import argparse
import glob
import lmfit
import os
import powerlaw
import random
import re
import string
import struct
import subprocess
import sys
import time
#
import networkx as nx
import numpy as np
import pandas as pd
import pickle as pk
import matplotlib.pyplot as plt
#
from collections import Counter, defaultdict, deque
from collections.abc import Iterable
from cycler import cycler
from dotenv import load_dotenv
from fractions import Fraction
from itertools import product
from networkx.classes.graph import Graph
from numpy.typing import NDArray
from networkx.drawing.layout import rescale_layout
from pathlib import Path
from os import chdir, getcwd
from os.path import join as pth_join
from scipy.cluster.hierarchy import fcluster, dendrogram, linkage
from scipy.interpolate import griddata, pchip
from scipy.linalg import expm, fractional_matrix_power
from scipy.linalg import eigvalsh as seigvalsh
from scipy.ndimage import gaussian_filter1d, gaussian_filter
from scipy.signal import argrelextrema, medfilt
from scipy.sparse import csr_array, spdiags, coo_matrix
from scipy.sparse import identity as scsp_identity
from scipy.sparse.linalg import eigsh as scsp_eigsh
from scipy.spatial.distance import squareform
from tqdm import tqdm
from typing import Any, Optional, Union, List, Tuple, Dict, Set, Type, \
    Sequence, Optional, Callable
