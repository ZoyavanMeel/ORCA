"""Module for plotting disparity curves."""

from typing import Optional
import numpy as np
import matplotlib.pyplot as plt


def plot_Z_curve_3D(x: np.ndarray, y: np.ndarray, z: np.ndarray, name: Optional[str] = None) -> None:
    """
    3D-plot function with name as title. `x`, `y`, and `z` should be the three components of the Z-curve.
    
    -----------------------------------------------------------------------
    Example
    >>> orca = ORCA.from_accession('NC_000913', 'example@email.com')
    >>> orca.find_oriCs()
    >>> plot_curves_3D(orca.x, orca.y, orca.z, name=orca.accession)
    """
    fig = plt.figure(figsize=(7,7))
    ax = plt.axes(projection='3d')
    ax.set_xlabel('Purine vs. Pyrimidine', fontsize=10, labelpad=15)
    ax.set_ylabel('Amino vs. Keto', fontsize=10, labelpad=15)
    ax.set_zlabel('Weak vs. Strong H-bond', fontsize=10, labelpad=15)

    ax.plot3D(x, y, z, c='b', linewidth=0.8)
    if name is not None:
        ax.set_title(f'Z-Curve: {name}', fontsize=10, loc='center', pad=20)
    plt.show()


def plot_curves(curves: tuple[np.ndarray], labels: list[str], peaks: Optional[list[int]], name: Optional[str] = None) -> None:
    """
    Plots up to 6 different y-axes onto a single figure. Ideal for displaying multiple disparity curves in a single plot.
    If displaying more than 3 diferent axes at once, some manual adjustment of the subplot paramenters might be needed.
    - `curves` : list of lists with y-axis values.
    - `labels` : list of names of each curve in `curves`.
    - `peaks`  : optional, list with indeces to plot onto the `curves`.
    - `name`   : optional, used in plot title.

    -----------------------------------------------------------------------
    Example
    >>> orca = ORCA.from_accession('NC_000913', 'example@email.com')
    >>> orca.find_oriCs()
    >>> plot_curves(curves=[orca.x, orca.y, orca.gc], labels=['$x_n$', '$y_n$', '$GC_n$'], peaks=orca.oriC_middles, name=orca.accession)

    Alternatively call:

    >>> orca = ORCA.from_accession('NC_000913', 'example@email.com')
    >>> orca.find_oriCs()
    >>> orca.plot_oriC_curves()
    """
    x_len = str(len(curves[0]))
    if int(x_len[1]) <= 4:
        x_max = x_len[0] + str(int(x_len[1])+(5-int(x_len[1]))) + '0'*len(x_len[2:])
    else:
        x_max = str(int(x_len[0])+1) + '0'*len(x_len[1:])
    thing = max(int(x_max)//1000, 1)
    xthing = thing * 100
    ything = thing*2

    # Colourblind friendly list of colours, I think.
    color_list = ['#e31a1c', '#1f77b4', '#33a02c', '#6a3d9a', '#ff7f00', '#b15928']

    fig = plt.figure(figsize=(8.5,4))
    fig.subplots_adjust(right=0.85 - (0.1 * max(len(curves)-2, 0)), bottom=0.25)
    base_ax = plt.axes()
    ax_list = [base_ax] + [base_ax.twinx() for i in range(len(curves) - 1)]

    offset = 1
    for axis in ax_list[1:]:
        axis.spines.right.set_position(("axes", offset))
        axis.yaxis.get_offset_text().set_x(offset)
        offset += 0.2

    good = False
    while not good:
        yticks_len = 0
        for i, ax in enumerate(curves):
            ubound, lbound = ything * round(max(ax)/ything), ything * round(min(ax)/ything)
            upper = ubound if max(ax) <= ubound else ubound + ything
            lower = lbound if min(ax) >= lbound else lbound - ything
            if len([x for x in range(lower, upper+ything, ything)]) > yticks_len:
                yticks_len = len([x for x in range(lower, upper+ything, ything)])
        if yticks_len < 6:
            yticks_len *= 2
            ything = ything // 2
        else:
            good = True
            break

    handle_list = []
    for i, ax in enumerate(curves):
        peaks_y = np.asarray([ax[j] for j in peaks]) # y refers to the y-axis coordinates, not the y-curve
        ax_list[i].plot(range(len(ax)), ax, color=color_list[i], zorder=2, label=labels[i])
        ax_list[i].scatter(peaks, peaks_y, marker='o', c='k', zorder=3, label='$\it{oriC}$')
        ax_list[i].tick_params(axis ='y', colors=color_list[i])
        ax_list[i].ticklabel_format(axis='y', style='sci', scilimits=(3, 3), useMathText=True)

        lbound = ything * round(min(ax)/ything)
        lower = lbound if min(ax) >= lbound else lbound - ything
        yticks = [lower + ything*j for j in range(yticks_len)]
        ax_list[i].set_yticks(yticks)
        ax_list[i].set_ylim(min(yticks), max(yticks))

        h, _ = ax_list[i].get_legend_handles_labels()
        handle_list.extend(h[:-1])
        oriC_handle = h[-1]
    handle_list.append(oriC_handle)

    if name is not None:
        base_ax.set_title(f'2D Z-Curve: {name}', fontsize=10,loc='center', pad=20)

    ubound= xthing * round(len(ax)/xthing)
    upper = ubound if len(ax) <= ubound else ubound + xthing
    xticks = [x for x in range(0, upper+xthing, xthing)]
    base_ax.set_xticks(xticks)

    base_ax.ticklabel_format(axis='x', style='sci', scilimits=(3,3), useMathText=True)
    base_ax.set_xlabel('Sequence length (bp)')
    base_ax.set_xlim(min(xticks), max(xticks))
    base_ax.grid(True, which='major', zorder=1)

    if len(peaks) != 0:
        l = labels + ['Prediction']
        n = len(curves)+1
    else:
        l = labels
        n = len(curves)

    plt.legend(
        handles=handle_list,
        labels=l,
        bbox_to_anchor=(0.12, -0.35, 0.75, .102),
        loc='center',
        ncol=n,
        mode="expand",
        borderaxespad=0.
    )
    plt.show()


def plot_skew(skewArray: np.ndarray, peaks: list[int], name: str) -> None:
    """Plots single skew diagram and its peaks"""

    fig = plt.figure()
    ax1 = plt.axes()

    peaks_y = np.asarray([skewArray[i] for i in peaks])

    ax1.set_title(f'GC-skew: {name}', fontsize=10,loc='center', pad=20)
    ax1.plot(range(len(skewArray)), skewArray, 'r', zorder=1)
    ax1.scatter(peaks, peaks_y, marker='X', c='k', zorder=2)
    plt.show()
