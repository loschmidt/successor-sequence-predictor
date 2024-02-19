__author__ = "Pavel Kohout <xkohou15@vutbr.cz>"
__date__ = "2024/02/16"
__description__ = "Keep allocation for all scripts in this study"

from pathlib import Path
from typing import Dict, List, Union

# General types
LOG = str

# Predictor
AA = str
FEATURE_VAL = float
AA_FEATURE_MAP = Dict[AA, FEATURE_VAL]

AA_INDICES = Dict[str, AA_FEATURE_MAP]

TREE_PATHS = List[Union[str, Path]]
LINEAGE_SEQUENCES = List[List[str]]  # List of lists with symbols for one MSA column from root to extant sequence
