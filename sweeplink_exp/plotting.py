import matplotlib.pyplot as plt
import numpy as np
from . import config, results, metrics

def draw_master_plot(char_name, threshold, save_file, n_total):
    sel_list = config.SIM_MAIN_LIST
    char_config = config.get_char_config(char_name)
    char_values = char_config["values"]

    # 1. Load data
    aggregated_data = {}
    for char_val in char_values:
        aggregated_data[char_val] = results.get_data_per_s_accross_repeats(
            char_name=char_name,
            char_val=char_val,
            sel_list=config.SIM_MAIN_LIST_WITH_0,
            n_total=n_total,
        )

    # 2. Setup figure
    fig, axes = plt.subplots(3, 1, figsize=(6, 6), sharex=True)
    x_centers = np.arange(len(char_values))

    # Panel A: draw TPR and FPR
    ax1 = axes[0]
    # get TPR, FPR and metrics
    metrics_data = {}
    TPR_y_data_per_sel = {sel: []  for sel in sel_list}
    FPR_y_data = []
    for char_val in char_values:
        TPR_dict, FPR, metrics_per_val = metrics.get_TPR_per_sel_and_FPR(
            parsed_data=aggregated_data[char_val],
            sel_list=[0.0] + sel_list,
            threshold=config.THRESHOLDS["sweeplink"],
            score_is_pval=False,
            return_metrics=True
        )
        for sel in sel_list:
            TPR_y_data_per_sel[sel].append(TPR_dict[sel])
        FPR_y_data.append(FPR)
        metrics_data[char_val] = metrics_per_val
    # draw
    for i, sel in enumerate(sel_list):
        ax1.plot(x_centers, TPR_y_data_per_sel[sel], marker='o', color=config.SEL_COLORS[sel], 
                 linewidth=2.5, markersize=8, label=rf"TPR, $s = {sel}$", zorder=config.SEL_ZORDER[sel])

    # Plot FPR Line & Shaded Area
    ax1.fill_between(x_centers, 0, FPR_y_data, color='gray', alpha=0.2, zorder=1)
    ax1.plot(x_centers, FPR_y_data, color='black', linestyle='--', linewidth=2.5, label='FPR', zorder=1)

    ax1.set_ylim(-0.05, 1.05)
    ax1.set_ylabel("Power (TPR or FPR)")
    ax1.set_title(f"(A) Detection Power and FPR {char_config['title']}", loc='left')
    ax1.legend(loc="lower right", prop={"size": 8})

    # Panel B: relative error to predict selection
    ax2 = axes[1]

    box_width = 0.2
    offsets = {0.01: -box_width, 0.02: 0, 0.05: box_width}

    for sel in sel_list:
        pred_data = {val: metrics.get_error(metrics_data[val][sel]["raw_sig_preds"], sel) if len(metrics_data[val][sel]["raw_sig_preds"]) > 2 else [] for val in char_values}

        # Calculate exact x positions for this specific selection coefficient
        positions = x_centers + offsets[sel]

        # Draw the boxplot (manage_ticks=False prevents it from messing up our X-axis)
        bplot = ax2.boxplot(
            [pred_data[val] for val in char_values],
            positions=positions,
            widths=box_width*0.8,
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
    ax2.axhline(0, color='black', linestyle='--', linewidth=1, zorder=0)

    ax2.set_ylabel(r"Relative Error of $\hat{s}$")
    ax2.set_title(r"(B) Estimation Accuracy of $\hat{s}$ " + f"{char_config['title']}", loc='left')

    # Panel C: violin plots of population size
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

    ax2.set_ylim(-1, 1)

    # True Ne reference line
    ax3.axhline(4.0, color='black', linestyle='--', linewidth=1, label='True $N_e$', zorder=0)


    ticks = [3, 5, 7, 9] # The log10 values you want to show
    ax3.set_yticks(ticks)
    ax3.set_yticklabels([f"$10^{i}$" for i in ticks])
    ax3.set_ylim(2, 10)

    ax3.set_ylabel("Posterior $N_e$")
    ax3.set_title(rf"(C) Posterior Uncertainty of $N_e$ {char_config['title']}", loc='left')
    ax3.legend(loc="upper right")

    ax3.set_xticks(x_centers)
    ax3.set_xticklabels(char_values)
    ax3.set_xlabel(char_config['x_label'])

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
