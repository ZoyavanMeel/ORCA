"""Script to visualise the dicision boundary"""

import os, sys
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import auc

MODEL_PATH = "data/output/machine_learning/24k_set_model.pkl.gz"
DATA_PATH  = "data/output/machine_learning/labels.csv"
PLOT_PATH  = "Machine_learning/DoriC_pairplot_06.png"

RANDOM_STATE = 42

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from _4_validate_train_model import load_data_labels_from_df, load_data_labels_from_path


def reduce_dimensions_and_plot_them():
    X, y = load_data_labels_from_path(DATA_PATH)
    start = time.perf_counter()
    X_norm = StandardScaler().fit_transform(X)
    print("Done with scaling: ", time.perf_counter() - start)
    start = time.perf_counter()
    X_2D = TSNE(random_state=RANDOM_STATE).fit_transform(X_norm)
    print("Done with reducing:", time.perf_counter() - start)
    start = time.perf_counter()

    data = np.concatenate([X_2D, y.to_numpy().reshape((-1, 1))], axis=1)
    data = pd.DataFrame(data, columns=['TSNE_x', 'TSNE_y', 'correct'])
    data['correct'].astype('int')

    # data.to_csv("Machine_learning/tsne.csv", index=False)
    # data = pd.read_csv('Machine_learning/tsne.csv')
    # print("Done with saving:", time.perf_counter() - start)

    start = time.perf_counter()
    plt.scatter(x=data[data['correct'] == 1]['TSNE_x'], y=data[data['correct'] == 1]['TSNE_y'], color='b', marker='.', s=0.1)
    sns.kdeplot(x=data[data['correct'] == 1]['TSNE_x'], y=data[data['correct'] == 1]['TSNE_y'], fill=True, alpha=0.7, cmap='Blues')
    print("Done with kde 1:", time.perf_counter() - start)
    start = time.perf_counter()
    plt.scatter(x=data[data['correct'] == 0]['TSNE_x'], y=data[data['correct'] == 0]['TSNE_y'], color='r', marker='.', s=0.1)
    sns.kdeplot(x=data[data['correct'] == 0]['TSNE_x'], y=data[data['correct'] == 0]['TSNE_y'], fill=True, alpha=0.5, cmap='Reds')
    print("Done with kde 2:", time.perf_counter() - start)
    start = time.perf_counter()
    plt.show()


def infer_cmap(color):
    """Apply the correct cmap to the chosen color; expand as needed."""
    if color == (0.413, 0.571, 1):
        return 'Blues'
    elif color == (1, 0.418, 0.425):
        return 'Reds'

def kde_hue(x, y, **kws):
    """Custom plot function to apply separate cmap coloring to each class in the classification problem."""
    ax = plt.gca()
    cmap = infer_cmap(kws['color'])
    sns.kdeplot(x=x, y=y, ax=ax, fill=True, cmap=cmap, **kws)
    return ax


def plot_pairgrid() -> None:
    """Plot every feature vs. every other feature in a pairplot"""

    start = time.perf_counter()
    print("Loading data............", end="")
    data = pd.read_csv(DATA_PATH)[['Z_score', 'G_score', 'D_score', 'total_pot', 'correct']]
    print(f"Done! ({time.perf_counter() - start:.3f})")
    start = time.perf_counter()

    color_dict = {
        True  : (0.413, 0.571, 1),
        False : (1, 0.418, 0.425)
    }

    print("Initialising PairGrid...", end="")
    g = sns.PairGrid(data, hue='correct', palette=color_dict, despine=False)
    print(f"Done! ({time.perf_counter() - start:.3f})")
    start = time.perf_counter()
    
    print("Plotting lower half.....", end="")
    g = g.map_lower(kde_hue, alpha=0.6)
    print(f"Done! ({time.perf_counter() - start:.3f})")
    start = time.perf_counter()

    print("Plotting diagonal.......", end="")
    g = g.map_diag(sns.kdeplot, fill=True)
    print(f"Done! ({time.perf_counter() - start:.3f})")
    start = time.perf_counter()

    print("Plotting upper half.....", end="")
    g = g.map_upper(plt.scatter, marker=".", s=0.2)
    print(f"Done! ({time.perf_counter() - start:.3f})")
    start = time.perf_counter()

    print("Final touches...........", end="")
    g = g.add_legend()

    for i, ax in enumerate(g.axes.flat):
        # set xticks
        if (i+1)%4 == 0:
            ax.set_xticks([i for i in range(-2,13, 2)])
            ax.set_xlim(-2, 12)
        else:
            ax.set_xticks([i/10-0.2 for i in range(0,15, 2)])
            ax.set_xlim(-0.2, 1.2)
        
        # set yticks
        if i not in [12, 13, 14, 15]:
            ax.set_yticks([i/10-0.2 for i in range(0,15, 2)])
            ax.set_ylim(-0.2, 1.2)
        else:
            ax.set_yticks([i for i in range(-2,13, 2)])
            ax.set_ylim(-2, 12)
    print(f"Done! ({time.perf_counter() - start:.3f})")
    start = time.perf_counter()
    
    print("Saving figure...........", end="")
    plt.savefig(PLOT_PATH, dpi=1200)
    print(f"Done! ({time.perf_counter() - start:.3f})")
    plt.show()


def precision_recall_plot() -> None:

    df = pd.read_csv("data/output/precision_recall/merged_precision_recall.csv").sort_values(by="threshold")

    fig, ax = plt.subplots()
    print("exp", auc(df["recall_exp"], df["precision_exp"]))
    print("30", auc(df["recall_30"], df["precision_30"]))
    ax.scatter(x=df["recall_exp"] * 100, y=df["precision_exp"] * 100, marker=".", c="r", label="experimental set")
    ax.scatter(x=df["recall_30"] * 100, y=df["precision_30"] * 100, marker=".", c="b", label="30% set")

    ax.set_axisbelow(True)
    ax.set(
        xlim=(0., 101), xticks=np.arange(0., 101, 10), xlabel="Recall (%)",
        ylim=(0., 101), yticks=np.arange(0., 101, 10), ylabel="Precision (%)"
    )

    ax.grid(True, which='major', color='k', linestyle='-')
    ax.legend(loc="lower left")
    # ax.grid(True, which='minor', color='grey', linestyle='--')
    # ax.grid(True, "both", zorder=0)

    # for i, txt in df["threshold"][:-1].items():
    #     ax.annotate(txt, (df[f"recall_{dataset}"][i], df[f"precision_{dataset}"][i]))

    plt.minorticks_on()
    plt.show()



if __name__ == "__main__":
    # plot_pairgrid()
    # reduce_dimensions_and_plot_them()
    precision_recall_plot()