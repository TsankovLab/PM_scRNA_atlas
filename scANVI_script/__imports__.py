""" -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : __imports__.py

* Description : All imports

* Creation Date : 02-20-2024

* Last Modified : Fri 08 Mar 2024 03:53:24 PM EST

* Created By : Atharva Bhagwat

_._._._._._._._._._._._._._._._._._._._._. """

import os
import scvi
import torch

import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

from pathlib import Path
from scvi.model.utils import mde
from scanpy._settings import ScanpyConfig

plt.rcParams['figure.dpi'] = 300
scvi.settings.seed = 42

USE_GPU = True if torch.cuda.is_available() else False
DEVICE = torch.device("cuda:0" if USE_GPU else "cpu")
