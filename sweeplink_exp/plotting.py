import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib.ticker import PercentFormatter
import numpy as np
from . import config, results, metrics

WIDTH = 0.6
HLINE_COLOR = "grey"
HLINE_WIDTH = 0.5
HLINE_ALPHA = 0.8


def draw_hline(ax, x_level):
    ax.axhline(x_level, color=HLINE_COLOR, linestyle='--', linewidth=HLINE_WIDTH, zorder=0, alpha=HLINE_ALPHA)

def get_offsets(n_bars, width=0.6):
    box_width = WIDTH / n_bars
    offsets = [(i - (n_bars - 1) / 2.0) * box_width for i in range(n_bars)]
    return offsets

def plot_power(ax, x_centers, sel_list, metrics_data, char_name, bar_plot=False):
    char_config = config.get_char_config(char_name)
    char_values = char_config["values"]

    TPR_y_data_per_sel = {sel: []  for sel in sel_list}
    FPR_y_data = []
    for char_val in char_values:
        tool_name = config.get_tool_name(char_name, char_val)
        TPR_dict, FPR = metrics.get_TPR_per_sel_and_FPR(
            metrics_data=metrics_data[char_val],
            sel_list=[0.0] + sel_list,
        )
        for sel in sel_list:
            if sel != 0.0:
                TPR_y_data_per_sel[sel].append(TPR_dict[sel])
        FPR_y_data.append(FPR)
    # draw
    if not bar_plot:
        for i, sel in enumerate(sel_list):
            if sel != 0:
                ax.plot(x_centers, TPR_y_data_per_sel[sel], marker='o', color=config.SEL_COLORS[sel],
                         linewidth=2.5, markersize=8, label=rf"$s = {sel}$", zorder=config.SEL_ZORDER[sel])
        # Plot FPR Line & Shaded Area
        ax.fill_between(x_centers, 0, FPR_y_data, color=config.SEL_COLORS[0.0], alpha=0.2, zorder=1)
        ax.plot(x_centers, FPR_y_data, color=config.SEL_COLORS[0.0], linestyle='--', linewidth=2.5, zorder=1)
        # Add the label
        ax.annotate('FPR',
            xy=(x_centers[-1], FPR_y_data[-1]),          # Anchor exactly to the last data point
            xytext=(0, 1),              # Offset text vertically by 1 point
            textcoords='offset points', # Tells matplotlib to use point offsets
            color='black',
            fontsize=8,
            va='bottom',                # Sits above the offset point
            ha='right')                 # Aligns it nicely to the right edge
        ax.set_ylim(-0.05, 1.05)
        ax.set_ylabel("Detection Rate (%)")
    else:
        offsets = get_offsets(n_bars=len(sel_list))
        min_shift = 0.02
        for sel_ind, sel in enumerate(sel_list):
            if sel == 0.0:
                y_data = FPR_y_data
            else:
                y_data = TPR_y_data_per_sel[sel]
            y_data = [_y + min_shift if (_y < min_shift and _y != 0) else _y for _y in y_data]
            ax.bar(
                x=x_centers + offsets[sel_ind],
                height=y_data,
                width=WIDTH/len(sel_list)*0.8,
                label=rf"$s = {sel}$",
                color=config.SEL_COLORS[sel],
                edgecolor='black',
                linewidth=1)
        fpr_container = ax.containers[0]
        custom_labels =[f"{(bar.get_height()-min_shift) * 100:.2f}%    " for bar in fpr_container]
        ax.bar_label(
            fpr_container,
            labels=custom_labels,  # Use our custom percentage strings
            padding=1,             # Adds a 3-point gap above the bar so it doesn't overlap
            fontsize=5,            # Make it easy to read
            color='black',
        )
        ax.set_ylabel("Detection Rate (%)")
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    for i in range(11):
        draw_hline(ax, x_level=i*0.1)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1.0, symbol=''))



def plot_boxplots(ax, x_centers, sel_list, metrics_data, char_values, relative_error=True):
    if relative_error:
        sel_list = [sel for sel in sel_list if sel != 0.0]
    offsets = get_offsets(n_bars=len(sel_list))

    for sel_ind, sel in enumerate(sel_list):
        pred_data = {}
        for val in char_values:
            if len(metrics_data[val][sel]["raw_sig_preds"]) > 2:
                if relative_error:
                    pred_data[val] = metrics.get_error(metrics_data[val][sel]["raw_sig_preds"], sel)
                else:
                    pred_data[val] = metrics_data[val][sel]["raw_sig_preds"]
            else:
                pred_data[val] = []

        # Calculate exact x positions for this specific selection coefficient
        positions = x_centers + offsets[sel_ind]

        # Draw the boxplot (manage_ticks=False prevents it from messing up our X-axis)
        bplot = ax.boxplot(
            [pred_data[val] for val in char_values],
            positions=positions,
            widths=WIDTH/len(sel_list)*0.8,
            patch_artist=True,
            manage_ticks=False,
            showfliers=False
        )

        # Color the boxes and medians
        for patch in bplot['boxes']:
            patch.set_facecolor(config.SEL_COLORS[sel])
            patch.set_edgecolor('black')
            patch.set_linewidth(1)
        for median in bplot['medians']:
            median.set_color('black')
            median.set_linewidth(1.5)

    # Perfect estimation reference line
    if relative_error:
        ax.axhline(0, color='black', linestyle='--', linewidth=1, zorder=0)
        ax.set_ylabel(r"Relative Error")
        ax.set_ylim(-1, 1)
        draw_hline(ax, x_level=0)
        for i in range(4):
            draw_hline(ax, i*0.25)
            draw_hline(ax, -i*0.25)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
    else:
        for sel in sel_list:
            ax.axhline(sel, color=config.SEL_COLORS[sel], linestyle='--', linewidth=1, zorder=0)
            if sel != 0:
                draw_hline(ax, x_level=-sel)
            draw_hline(ax, x_level=0.06)
            draw_hline(ax, x_level=-0.06)
        ax.set_ylabel(r"Estimated Selection")
        ax.set_ylim(-0.065, 0.065)

        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.02))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.01))

def draw_master_plot(char_name, save_file, n_total):
    sel_list = config.get_sel_list(char_name)
    char_config = config.get_char_config(char_name)
    char_values = char_config["values"]

    # -1. Check if we have comparison of tools - we have to draw slightly different plot
    if char_name == "comparison":
        A_bar_plot = True
        B_relative_error = False
        C_draw = False
        fig_width = 8
    else:
        A_bar_plot = False
        B_relative_error = True
        C_draw = True
        fig_width = 6

    # 1. Load data
    metrics_data = {}
    for char_val in char_values:
        parsed_data = results.get_data_per_s_accross_repeats(
            char_name=char_name,
            char_val=char_val,
            n_total=n_total,
        )
        tool_name = config.get_tool_name(char_name, char_val)
        metrics_data[char_val] = metrics.calculate_metrics(
            parsed_data=parsed_data,
            sel_list=[0.0] + sel_list,
            threshold=config.THRESHOLDS[tool_name],
            score_is_pval=config.IS_PVAL[tool_name]
        )

    # 2. Setup figure
    fig, axes = plt.subplots(3 - int(not C_draw), 1, figsize=(fig_width, 2 * (3 - int(not C_draw))), sharex=True)
    x_centers = np.arange(len(char_values))

    # Panel A: draw TPR and FPR
    ax1 = axes[0]
    plot_power(ax1, x_centers, [0.0] + sel_list, metrics_data, char_name, bar_plot=A_bar_plot)

    ax1.set_title(f"(A) Sweep Detection Performance {char_config['title']}", loc='left', size=11)
    if char_name != "comparison":
        ax1.legend(loc="center right", prop={"size": 8})
    else:
        ax1.legend(
            loc="upper left",
            bbox_to_anchor=(1.05, 1), # Moves legend outside the plot
            borderaxespad=0.
        )

    # Panel B: relative error to predict selection
    ax2 = axes[1]
    plot_boxplots(ax2, x_centers, [0.0] + sel_list, metrics_data, char_values, relative_error=B_relative_error)
    ax2.set_title(r"(B) Accuracy of Detected Sweeps " + f"{char_config['title']}", loc='left', size=11)

    # Panel C: violin plots of population size
    if C_draw:
        ax3 = axes[2]

        posterior_samples = {val: results.get_trace_data_accross_repeats(char_name, val, n_total=100) for val in char_values}

        # Draw the violins exactly at the x_centers
        vplot = ax3.violinplot([posterior_samples[val]["log10N_pop0"] for val in char_values], positions=x_centers, showmedians=True, showextrema=True)

        # Format the violins (Neutral light gray)
        for pc in vplot['bodies']:
            pc.set_facecolor('lightgray')
            pc.set_edgecolor('black')
            pc.set_alpha(1) # Solid, not transparent

        # Format the median and range lines inside the violin
        for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
            vp = vplot[partname]
            vp.set_edgecolor('black')
            vp.set_linewidth(1)

        # True Ne reference line
        ax3.axhline(4.0, color='black', linestyle='--', linewidth=1, label='True $N_e$', zorder=0)


        ticks = [3, 5, 7, 9] # The log10 values you want to show
        ax3.set_yticks(ticks)
        ax3.set_yticklabels([f"$10^{i}$" for i in ticks])
        ax3.set_ylim(2, 10)

        ax3.set_ylabel("Posterior $N_e$")
        ax3.set_title(rf"(C) Posterior Uncertainty of $N_e$ {char_config['title']}", loc='left', size=11)
        ax3.legend(loc="upper right")

    axes[-1].set_xticks(x_centers)
    if char_name != "comparison":
        axes[-1].set_xticklabels(char_values)
        axes[-1].set_xlabel(char_config['x_label'])
    else:
        axes[-1].set_xticklabels([config.TOOL_LABELS[char_val] for char_val in char_values])

    default_value = char_config['default_value']
    if default_value is not None:
        default_idx = char_values.index(default_value)
        for ax in axes: # Apply to Panel A, B, and C
            ax.axvspan(default_idx - 0.5, default_idx + 0.5,
                       color='#E6F3FF', alpha=0.8, zorder=-1)

    # Ensure panels don't overlap
    plt.tight_layout()

    if save_file is not None:
        plt.savefig(save_file)
    else:
        plt.show()


def draw_MCC_plot(char_name, save_file, n_total):
    char_values = config.get_char_config(char_name)["values"]
    sel_list = config.get_sel_list(char_name)

    # 1. Load data
    aggregated_data = {}
    for char_val in char_values:
        aggregated_data[char_val] = results.get_data_per_s_accross_repeats(
        char_name=char_name,
            char_val=char_val,
            n_total=n_total,
        )

    # 2. Create grid and figure
    fig, ax = plt.subplots(figsize=(4, 2.7))
    ax.set_ylim([-0.05, 1.05])
    limit = 8
    grid = np.concatenate([[0], np.geomspace(10**(-limit), 1.0, 500)])

    x_axis = grid

    best_results = {}
    target_fpr = None
    target_fpr_results = {}

    for i, char_val in enumerate(char_values):
        tool_name = config.get_tool_name(char_name, char_val)
        y_data = []
        fpr_data = []
        for thr in grid:
            thr_upd = thr
            if not config.IS_PVAL[tool_name]:
                thr_upd = 1 - thr
            mcc_score, fpr = metrics.get_MCC(
                parsed_data=aggregated_data[char_val],
                sel_list=[0.0] + sel_list,
                threshold=thr_upd,
                score_is_pval=config.IS_PVAL[tool_name],
                return_FPR=True,
            )
            y_data.append(mcc_score)
            fpr_data.append(fpr)
        if char_name == "comparison":
            color = config.TOOL_COLORS[tool_name]
            zorder = config.TOOL_ZORDER[tool_name]
            label = config.TOOL_LABELS[tool_name]
        else:
            color = config.TAB_CMAP(i)
            zorder = i
            label = f"{char_name}: {char_val}"
        plt.plot(x_axis, y_data, color=color, zorder=zorder, label=label)

        ind = np.argmax(y_data)
        best_results[tool_name] = [x_axis[ind], y_data[ind]]
        if target_fpr is None:
            target_fpr = fpr_data[ind]
        fpr_ind = np.argmin(np.abs(np.array(fpr_data) - target_fpr))
        target_fpr_results[tool_name] = [x_axis[fpr_ind], fpr_data[fpr_ind]]
        
        plt.plot(x_axis[ind], y_data[ind], color=color, zorder=zorder, marker="o", markersize=5)

        if config.IS_PVAL[tool_name]:
            p_string = f"pval={x_axis[ind]: .1e}"
        else:
            p_string = r"P=" + f"{1 - x_axis[ind]: .2f}"
        ax.annotate(f"MCC={y_data[ind]:.2f}\n{p_string}",
                    xy=(x_axis[ind], y_data[ind]),
                    xytext=(-4, 2),
                    textcoords='offset points',
                    fontsize=6,
                    color=color,
                    fontweight='bold',
                    ha='right',       # <--- ALIGN TEXT TO THE RIGHT
                    va='bottom',      # Vertical alignment
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))

    for tool in best_results:
        print(tool, f"best MCC={best_results[tool][1]} at threshold={best_results[tool][0]}")
    
    for tool in target_fpr_results:
        print(tool, f"FPR={target_fpr_results[tool][1]} at threshold={target_fpr_results[tool][0]}")

    plt.xscale('symlog', linthresh=10 ** (-limit))
    ticks = [1.0, 1e-1, 1e-3, 1e-5, 1e-7, 0.0]
    labels = ["1", r"$10^{-1}$", r"$10^{-3}$", r"$10^{-5}$", r"$10^{-7}$", "0"]
    plt.xticks(ticks, labels=labels)

    plt.xlim(1.05, 0.0)
    ax.set_xlabel(r"Significance Threshold (p-value or $1−P$)")
    ax.set_ylabel("Detection Power (MCC)")

    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=False)

    plt.tight_layout()

    if save_file is not None:
        plt.savefig(save_file)
    else:
        plt.show()
