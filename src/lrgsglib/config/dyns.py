from random import sample
import subprocess
from tqdm import tqdm

from ..nx_patches.funcs import *
from ..nx_patches.objects import *
from .utils import *
from .ds.BinDynSys import BinDynSys
from .ds.IsingDynamics import IsingDynamics, IsingDynamics_DEV
from .ds.SignedRW import SignedRW
from .ds.VoterModel import VoterModel
