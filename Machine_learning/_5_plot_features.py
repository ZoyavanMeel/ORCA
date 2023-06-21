"""Script to visualise the dicision boundary"""

import os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import joblib
from sklearn.decomposition import PCA
from sklearn.inspection import DecisionBoundaryDisplay

MODEL_PATH = "data/output/machine_learning/24k_set_model.pkl"
DATA_PATH  = "data/output/machine_learning/labels.csv"

RANDOM_STATE = 42

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from _4_validate_train_model import load_data_labels


def decompose_2D(path: str) -> tuple[pd.DataFrame, pd.Series]:
    X, y = load_data_labels(path)
    X_2D = PCA(n_components=2, random_state=RANDOM_STATE).fit_transform(X)
    return X_2D, y


def plot_xy_line(X: pd.DataFrame, y: pd.Series, model) -> None:
    disp = DecisionBoundaryDisplay.from_estimator(
        estimator=model,
        X=X,
        response_method="predict",
        alpha=0.5
    )
    disp.ax_.scatter(X[:, 0], X[:, 1], c=y, edgecolor="k", marker='.')
    plt.show()


def main() -> None:
    model = joblib.load(MODEL_PATH)
    X, y = decompose_2D(DATA_PATH)
    plot_xy_line(X, y, model)



if __name__ == "__main__":
    main()