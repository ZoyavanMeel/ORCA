"""Script to visualise the dicision boundary"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.inspection import DecisionBoundaryDisplay

from BioFile import comp_path

EXP_DATA = comp_path("data\output\machine_learning\labels_against_NCBI.csv")
DORIC_DATA = comp_path("data\output\machine_learning\labels_against_DoriC.csv")


def plot_xy_line(data: pd.DataFrame, line: np.ndarray) -> None:
    """
    `data`: 2-column DataFrame with the coordinates of the samples
    `line`: 2-column ndarray for plotting the decision boundary
    """
    DecisionBoundaryDisplay.from_estimator()


def main() -> None:
    ...


if __name__ == "__main__":
    main()