import matplotlib.pyplot as plt
import numpy as np

def plot_Z_curve_3D(Z_curve: tuple, name: str = None):
    """
    3D-plot function with name as title
    
    -----------------------------------------------------------------------
    Example
    >>> properties = find_oriCs(accession='NC_000913', email='example@email.com')
    >>> plot_Z_curve_3D(Z_curve=properties['z_curve'], name='NC_000913')
    """
    fig = plt.figure(figsize=(7,7))
    ax = plt.axes(projection='3d')
    ax.set_xlabel('Purine vs. Pyrimidine', fontsize=10, labelpad=15)
    ax.set_ylabel('Amino vs. Keto', fontsize=10, labelpad=15)
    ax.set_zlabel('Weak vs. Strong H-bond', fontsize=10, labelpad=15)

    x, y, z = Z_curve
    ax.plot3D(x, y, z, c='b', linewidth=0.8)
    if name is not None:
        ax.set_title(f'Z-Curve: {name}', fontsize=10, loc='center', pad=20)
    plt.show()


def plot_Z_curve_2D(curves, peaks, labels, name=None):
    """
    Plots 2D Z-curve. Can display up to 4 y-axes in a single figure.
    - `curves` : list of lists with y-axis values
    - `peaks`  : list with indeces of peaks for arrays in y_val_list
    - `name`   : used in plot title

    -----------------------------------------------------------------------
    Example
    >>> properties = find_oriCs(accession='NC_000913', email='example@email.com')
    >>> x_comp, y_comp = properties['z_curve'][:2]
    >>> gc_skew = properties['gc_skew']
    >>> oriCs = properties['oriC_middles']
    >>> plot_Z_curve_2D(curves=[x_comp, y_comp, gc_skew], peaks=oriCs, labels=['$x_n$', '$y_n$', '$g_n$'], name='NC_000913')
    """
    x_len = str(len(curves[0]))
    if int(x_len[1]) <= 4:
        x_max = x_len[0] + str(int(x_len[1])+(5-int(x_len[1]))) + '0'*len(x_len[2:])
    else:
        x_max = str(int(x_len[0])+1) + '0'*len(x_len[1:])
    thing = int(x_max)//1000
    xthing = thing * 100
    ything = thing*2

    color_list = ['r', 'b', 'g', 'c']
    fig = plt.figure(figsize=(8,4))
    fig.subplots_adjust(right=0.75, bottom=0.25)
    base_ax = plt.axes()
    ax_list = [base_ax] + [base_ax.twinx() for i in range(len(curves) - 1)]

    offset = 1
    for axis in ax_list[1:]:
        axis.spines.right.set_position(("axes", offset))
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
        ax_list[i].ticklabel_format(axis='y', style='sci', useMathText=True)

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

    plt.legend(
        handles=handle_list,
        labels=labels + ['$\it{oriC}$'],
        bbox_to_anchor=(0.12, -0.35, 0.75, .102),
        loc='center',
        ncol=len(curves)+1,
        mode="expand",
        borderaxespad=0.
    )
    plt.show()


def plot_skew(skewArray, peaks, name):
    """Plots single skew diagram and its peaks"""

    fig = plt.figure()
    ax1 = plt.axes()

    peaks_y = np.asarray([skewArray[i] for i in peaks])

    ax1.set_title(f'GC-skew: {name}', fontsize=10,loc='center', pad=20)
    ax1.plot(range(len(skewArray)), skewArray, 'r', zorder=1)
    ax1.scatter(peaks, peaks_y, marker='X', c='k', zorder=2)
    plt.show()


def distance_histogram(db, log=False):
    plt.hist(db['Distance_bp'], bins=[x for x in range(0, 3600000//2, 1000)], log=log) # 3.6e6 = avg. len of bacterial chromosome
    plt.show()
    plt.hist(db['Distance_pc'], bins=[x for x in range(0, 50, 1)], log=log)
    plt.show()