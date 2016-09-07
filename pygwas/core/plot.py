import scipy as sp
import matplotlib
import matplotlib.pyplot as plt

def simple_log_qqplot(quantiles_list, png_file=None, pdf_file=None, quantile_labels=None, line_colors=None,
            max_val=5, title=None, text=None, plot_label=None, ax=None, **kwargs):
    storeFig = False
    if ax is None:
        f = plt.figure(figsize=(5.4, 5))
        ax = f.add_axes([0.1, 0.09, 0.88, 0.86])
        storeFig = True
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=2.0)
    num_dots = len(quantiles_list[0])
    exp_quantiles = sp.arange(1, num_dots + 1, dtype='single') / (num_dots + 1) * max_val
    for i, quantiles in enumerate(quantiles_list):
        if line_colors:
            c = line_colors[i]
        else:
            c = 'b'
        if quantile_labels:
            ax.plot(exp_quantiles, quantiles, label=quantile_labels[i], c=c, alpha=0.5, linewidth=2.2)
        else:
            ax.plot(exp_quantiles, quantiles, c=c, alpha=0.5, linewidth=2.2)
    ax.set_ylabel("Observed $-log_{10}(p$-value$)$")
    ax.set_xlabel("Expected $-log_{10}(p$-value$)$")
    if title:
        ax.title(title)
    max_x = max_val
    max_y = max(map(max, quantiles_list))
    ax.axis([-0.025 * max_x, 1.025 * max_x, -0.025 * max_y, 1.025 * max_y])
    if quantile_labels:
        fontProp = matplotlib.font_manager.FontProperties(size=10)
        ax.legend(loc=2, numpoints=2, handlelength=0.05, markerscale=1, prop=fontProp, borderaxespad=0.018)
    y_min, y_max = plt.ylim()
    if text:
        f.text(0.05 * max_val, y_max * 0.9, text)
    if plot_label:
        f.text(-0.138 * max_val, y_max * 1.01, plot_label, fontsize=14)
    if storeFig == False:
        return
    if png_file != None:
        f.savefig(png_file)
    if pdf_file != None:
        f.savefig(pdf_file, format='pdf')


def simple_qqplot(quantiles_list, png_file=None, pdf_file=None, quantile_labels=None, line_colors=None,
            title=None, text=None, ax=None, plot_label=None, **kwargs):
    storeFig = False
    if ax is None:
        f = plt.figure(figsize=(5.4, 5))
        ax = f.add_axes([0.11, 0.09, 0.87, 0.86])
        storeFig = True
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=2.0)
    num_dots = len(quantiles_list[0])
    exp_quantiles = sp.arange(1, num_dots + 1, dtype='single') / (num_dots + 1)
    for i, quantiles in enumerate(quantiles_list):
        if line_colors:
            c = line_colors[i]
        else:
            c = 'b'
        if quantile_labels:
            ax.plot(exp_quantiles, quantiles, label=quantile_labels[i], c=c, alpha=0.5, linewidth=2.2)
        else:
            ax.plot(exp_quantiles, quantiles, c=c, alpha=0.5, linewidth=2.2)
    ax.set_ylabel("Observed $p$-value")
    ax.set_xlabel("Expected $p$-value")
    if title:
        ax.title(title)
    ax.axis([-0.025, 1.025, -0.025, 1.025])
    if quantile_labels:
        fontProp = matplotlib.font_manager.FontProperties(size=10)
        ax.legend(loc=2, numpoints=2, handlelength=0.05, markerscale=1, prop=fontProp, borderaxespad=0.018)
    if text:
        f.text(0.05, 0.9, text)
    if plot_label:
        f.text(-0.151, 1.04, plot_label, fontsize=14)
    if storeFig == False:
        return
    if png_file != None:
        f.savefig(png_file)
    if pdf_file != None:
        f.savefig(pdf_file, format='pdf')


def plot_simple_qqplots(png_file_prefix, results, result_labels=None, line_colors=None,
            num_dots=1000, title=None, max_neg_log_val=5):
    """
    Plots both log QQ-plots and normal QQ plots.
    """
    qs = []
    log_qs = []
    for res in results:
        pvals = res.snp_results['scores'][:]
        qs.append(get_quantiles(pvals, num_dots))
        log_qs.append(get_log_quantiles(pvals, num_dots, max_neg_log_val))
    simple_qqplot(qs, png_file_prefix + '_qq.png', quantile_labels=result_labels,
                line_colors=line_colors, num_dots=num_dots, title=title)
    simple_log_qqplot(log_qs, png_file_prefix + '_log_qq.png', quantile_labels=result_labels,
                line_colors=line_colors, num_dots=num_dots, title=title, max_val=max_neg_log_val)


def plot_simple_qqplots_pvals(png_file_prefix, pvals_list, result_labels=None, line_colors=None,
            num_dots=1000, title=None, max_neg_log_val=5):
    """
    Plots both log QQ-plots and normal QQ plots.
    """
    qs = []
    log_qs = []
    for pvals in pvals_list:
        qs.append(get_quantiles(pvals, num_dots))
        log_qs.append(get_log_quantiles(pvals, num_dots, max_neg_log_val))
    simple_qqplot(qs, png_file_prefix + '_qq.png', quantile_labels=result_labels,
                line_colors=line_colors, num_dots=num_dots, title=title)
    simple_log_qqplot(log_qs, png_file_prefix + '_log_qq.png', quantile_labels=result_labels,
                line_colors=line_colors, num_dots=num_dots, title=title, max_val=max_neg_log_val)


