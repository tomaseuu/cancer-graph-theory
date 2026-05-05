import random
from pathlib import Path

import numpy as np


def ensure_output_dir(path):
    output_path = Path(path)
    directory = output_path.parent if output_path.suffix else output_path
    directory.mkdir(parents=True, exist_ok=True)


def set_random_seed(seed=42):
    random.seed(seed)
    np.random.seed(seed)
