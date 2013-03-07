import sys
sys.path.append('../')

import numpy as np
from src.datawarehouse import DataWarehouse as Dw
from src.patch import Patch
from src.material import Material
from src.materialmodel2d import MaterialModel
from src.shape2 import GIMP as Shape
from src.boundcond import BoundaryCondition as Bc
from src import geomutils
import src.mpm2d as mpm