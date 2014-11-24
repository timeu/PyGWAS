    def plot_histogram(self, title=None , pdf_file=None, png_file=None, x_label=None, p_her=None,
            p_her_pval=None):

        if title:
            plt.figure(figsize=(6, 5.4))
            plt.axes([0.13, 0.11, 0.85, 0.82])
        else:
            plt.figure(figsize=(6, 4.8))
            plt.axes([0.13, 0.11, 0.85, 0.86])
        if x_label:
            plt.xlabel(x_label)
        phen_vals = self.get_values

        minVal = min(phen_vals)
        maxVal = max(phen_vals)
        x_range = maxVal - minVal
        histRes = plt.hist(phen_vals, bins=round(8 + 2 * sp.log(len(phen_vals))), alpha=0.7)
        y_max = max(histRes[0])
        plt.axis([minVal - 0.035 * x_range, maxVal + 0.035 * x_range, -0.035 * y_max, 1.19 * y_max])
        num_phen_vals = len(phen_vals)
        shapiro_pval = sp.stats.shapiro(phen_vals)[1]
        if p_her:
            if p_her_pval != None:
                st = "Num. of values: %d, herit.: %0.4f, herit. -log(p): %0.4f, transf.: %s" % \
                    (num_phen_vals, p_her, -sp.log10(p_her_pval), str(self.phen_dict[pid]['transformation']))
                plt.text(maxVal - 0.95 * x_range, 1.1 * y_max, st, size="xx-small")
            else:
                st = "Number of values: %d,  Pseudo-heritability: %0.4f,  Transformation: %s" % \
                    (num_phen_vals, p_her, str(self.phen_dict[pid]['transformation']))
                plt.text(maxVal - 0.95 * x_range, 1.1 * y_max, st, size="xx-small")
        else:
            st = "Number of values: %d, Transformation: %s" % (num_phen_vals, str(self.phen_dict[pid]['transformation']))
            plt.text(maxVal - 0.9 * x_range, 1.1 * y_max, st, size="x-small")
        plt.text(maxVal - 0.85 * x_range, 1.02 * y_max, "Shapiro-Wilk normality $p$-value: %0.6f" % shapiro_pval , size="x-small")
        #print max(histRes[0])
        plt.ylabel("Frequency")
        if title:
            plt.title(title)
        if pdf_file:
            plt.savefig(pdf_file, format="pdf")
        if png_file:
            plt.savefig(png_file, format="png", dpi=300)
        elif not pdf_file or png_file:
            plt.show()
        plt.clf()




    def plot_marker_box_plot(self, pid, marker, m_accessions, m_position=None, m_chromosome=None, plot_file=None,
                plot_format='png', title=None, m_score=None):
        """
        Plots a box plot for the given binary marker and phenotype. 
        
        Assumes the marker is integer based.        
        Assumes the marker and the phenotype accessions are aligned.
        """
        phen_vals = self.get_values(pid)
        if len(m_accessions) != len(phen_vals):
            raise Exception

        nt_counts = sp.bincount(marker)
        if len(nt_counts) > 2:
            import warnings
            warnings.warn("More than 2 alleles, box-plot might be wrong?")

        allele_phen_val_dict = {}
        for nt in set(marker):
            allele_phen_val_dict[nt] = {'values':[], 'ecotypes':[]}

        for i, nt in enumerate(marker):
            allele_phen_val_dict[nt]['values'].append(phen_vals[i])
            if m_accessions:
                allele_phen_val_dict[nt]['ecotypes'].append(m_accessions[i])

        xs = []
        positions = []
        for nt in allele_phen_val_dict:
            positions.append(nt)
            xs.append(allele_phen_val_dict[nt]['values'])
        plt.figure()
        plt.boxplot(xs, positions=positions)
        min_val = min(phen_vals)
        max_val = max(phen_vals)
        val_range = max_val - min_val
        max_pos = max(positions)
        min_pos = min(positions)
        x_range = max_pos - min_pos
        plt.axis([min_pos - 0.5 * x_range, max_pos + 0.5 * x_range, min_val - val_range * 0.3, max_val + val_range * 0.3])
        plt.text(min_pos - 0.45 * x_range, min_val - 0.15 * val_range, "# of obs.: ", color='k')
        for i, (x, pos) in enumerate(it.izip(xs, positions)):
            plt.text(pos - 0.05, min_val - 0.15 * val_range, str(len(xs[i])), color='k')
        if m_score:
            plt.text(min_pos + 0.13 * x_range, max_val + 0.15 * val_range,
                '$-log_{10}$(p-value)/score: %0.2f' % m_score, color='k')
        if title:
            plt.title(title)
        elif m_chromosome and m_position:
            plt.title('%s : chromosome=%d, position=%d' % (self.get_name(pid), m_chromosome, m_position))
        if plot_file:
            plt.savefig(plot_file, format=plot_format)
        else:
            plt.show()
        plt.clf()


    def plot_marker_accessions_hist(self, pid, marker, m_accessions, plot_file=None, plot_format='png',
                m_position=None, m_chromosome=None, title=None, m_score=None):
        """
        A histogram displaying the phenotype values (ordered) on the y-axis, and the accession on the x-axis.
        """
        import matplotlib.cm as cm
        import matplotlib.colors as colors

        color_map = {}
        colors = ['r', 'm', 'b', 'g']
        proxy_rects = []
        labels = []
        for nt in set(marker):
            c = colors.pop()
            color_map[nt] = c
            proxy_rects.append(plt.Rectangle((0, 0), 1, 1, fc=c, alpha=0.6))
            labels.append("'%s' allele" % str(nt))


        phen_values = self.get_values(pid)
        ets = self.get_ecotypes(pid)
        l = zip(phen_values, ets, marker)
        l.sort(reverse=True)

        fig = plt.figure(figsize=(20, 10))
        ax = fig.add_axes([0.07, 0.15, 0.91, 0.82])
        x_range = len(l) - 0.2
        min_y = min(phen_values)
        max_y = max(phen_values)
        y_range = max_y - min_y

        for i, (phen_val, accession, nt) in enumerate(l):
            color = color_map[nt]
            rect = ax.bar(i, phen_val, 0.8, color=color, alpha=0.6)
        ax.axis([-x_range * 0.02, x_range * 1.02, min_y - 0.05 * y_range, max_y + 0.05 * y_range])

        ax.legend(proxy_rects, labels)
        ax.set_ylabel('Phenotype value')
        ax.set_xticks((sp.arange(len(ets)) + 0.4).tolist())
        ax.set_xticklabels(ets, rotation='vertical', fontsize='xx-small')
        ax.set_xlabel('Ecotype IDs')

        fig.savefig(plot_file, format=plot_format, dpi=300)


    def plot_phen_relatedness(self, k, k_accessions, plot_file_prefix, pids=None):
        import linear_models as lm
        import pylab
        import scipy as sp
        from scipy import linalg
        if not pids:
            pids = self.get_pids()
        self.convert_to_averages(pids)
        self.filter_ecotypes_2(k_accessions, pids)
        for pid in pids:
            ets = self.get_ecotypes(pid)
            vals = self.get_values(pid)
            k_m = lm.prepare_k(k, k_accessions, ets)
            c = sp.sum((sp.eye(len(k_m)) - (1.0 / len(k_m)) * sp.ones(k_m.shape)) * sp.array(k_m))
            k_scaled = (len(k) - 1) * k / c
            p_her = self.get_pseudo_heritability(pid, k_m)
            x_list = []
            y_list = []
            for i in range(len(ets)):
                for j in range(i):
                    x_list.append(k_m[i, j])
                    y_list.append(vals[i] - vals[j])
            ys = sp.array(y_list)
            ys = ys * ys
            xs = sp.array(x_list)
            phen_name = self.get_name(pid)
            phen_name = phen_name.replace('<i>', '')
            phen_name = phen_name.replace('</i>', '')
            phen_name = phen_name.replace('+', '_plus_')
            phen_name = phen_name.replace('/', '_div_')
            file_name = plot_file_prefix + '_%d_%s.png' % (pid, phen_name)
            pylab.figure()
            pylab.plot(xs, ys, 'k.', alpha=0.2)
            pylab.xlabel('Relatedness')
            pylab.ylabel('Squared phenotypic difference')
            #Plot regression line
            Y_mat = sp.mat(ys).T
            X_mat = sp.hstack((sp.mat(sp.ones(len(xs))).T, sp.mat(xs).T))
            (betas, residues, rank, s) = linalg.lstsq(X_mat, Y_mat)
            x_min, x_max = pylab.xlim()
            pylab.plot([x_min, x_max], [betas[0] + x_min * betas[1], betas[0] + x_max * betas[1]])
            corr = sp.corrcoef(xs, ys)[0, 1]
            y_min, y_max = pylab.ylim()
            x_range = x_max - x_min
            y_range = y_max - y_min
            pylab.axis([x_min - 0.025 * x_range, x_max + 0.025 * x_range,
                    y_min - 0.025 * y_range, y_max + 0.15 * y_range])
            pylab.text(x_min + 0.1 * x_range, y_max + 0.03 * y_range, 'Correlation: %0.4f' % (corr))
            pylab.text(x_min + 0.5 * x_range, y_max + 0.03 * y_range, 'Pseudo-heritability: %0.4f' % (p_her))
            pylab.savefig(file_name)
            del k_m
            del k_scaled










def plot_snp_map(self, chromosome, position, pdf_file=None, png_file=None, map_type='global',
            color_by=None, cmap=None, title='', eid=None):
        """
        Plot accessions on a map.
        
        'color_by' is by default set to be the phenotype values.
        """
        import phenotypeData as pd
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        #matplotlib.rcParams['backend'] = 'GTKAgg'
        if eid == None:
            eid = pd.get_ecotype_id_info_dict()
        lats = []
        lons = []
        acc_names = []
        for e in self.accessions:
            r = eid[int(e)]
            acc_names.append(r[0])
            try:
                latitude = float(r[2])
                longitude = float(r[3])
#                r = eid[str(e)]
#                latitude = float(r[5])
#                longitude = float(r[6])

            except Exception, err_str:
                print "Latitude and Longitude, not found?:", err_str
                print 'Placing them in the Atlantic.'
                latitude = 40
                longitude = -20

            lats.append(latitude)
            lons.append(longitude)

        from mpl_toolkits.basemap import Basemap
        import numpy as np
        from pylab import cm
        if map_type == "global2":
            plt.figure(figsize=(14, 12))
            m = Basemap(width=21.e6, height=21.e6, projection='gnom', lat_0=76, lon_0=15)
            m.drawparallels(np.arange(20, 90, 20))
            m.drawmeridians(np.arange(-180, 180, 30))
        elif map_type == 'global':

            plt.figure(figsize=(16, 4))
            plt.axes([0.02, 0.02, 0.96, 0.96])
            m = Basemap(projection='cyl', llcrnrlat=10, urcrnrlat=80,
                    llcrnrlon= -130, urcrnrlon=150, lat_ts=20, resolution='c')
            m.drawparallels(np.arange(20, 90, 20))
            m.drawmeridians(np.arange(-180, 180, 30))
        elif map_type == 'europe':

            plt.figure(figsize=(8, 6))
            plt.axes([0.02, 0.02, 0.96, 0.96])
            m = Basemap(projection='cyl', llcrnrlat=35, urcrnrlat=70,
                    llcrnrlon= -15, urcrnrlon=40, lat_ts=20, resolution='h')
            m.drawparallels(np.arange(30, 80, 10))
            m.drawmeridians(np.arange(-20, 100, 10))
            #m.bluemarble()
        elif map_type == 'sweden':

            plt.figure(figsize=(2.4, 4))
            plt.axes([0.02, 0.02, 0.96, 0.96])
            m = Basemap(projection='merc', llcrnrlat=55, urcrnrlat=67,
                    llcrnrlon=10, urcrnrlon=25, lat_ts=10, resolution='i')
            m.drawparallels(np.arange(45, 75, 5))
            m.drawmeridians(np.arange(5, 30, 5))
            #m.bluemarble()
        else:
            raise Exception("map_type is invalid")

        #m.drawmapboundary(fill_color='aqua')
        m.drawcoastlines()
        m.fillcontinents()
        m.drawcountries()
        #m.fillcontinents(color='green', lake_color='blue')

        xs = []
        ys = []
        for lon, lat in zip(lons, lats):
            x, y = m(*np.meshgrid([lon], [lat]))
            xs.append(float(x))
            ys.append(float(y))

        if not color_by:
            color_vals = self.get_snp_at(chromosome, position)
        else:
            color_vals = color_by
        assert len(color_vals) == len(self.accessions), "accessions and color_by_vals values don't match ! "
        if not cmap:
            num_colors = len(set(color_vals))
            if num_colors <= 10:
                cmap = cm.get_cmap('jet', num_colors)
            else:
                cmap = cm.get_cmap('jet')
        lws = [0] * len(xs)
        plt.scatter(xs, ys, s=10, linewidths=lws, c=color_vals, cmap=cmap, alpha=0.7, zorder=2)
        #plt.plot(xs, ys, 'o', color='r', alpha=0.5, zorder=2,)
        if title:
            plt.title(title)
        if pdf_file:
            plt.savefig(pdf_file, format="pdf", dpi=400)
        if png_file:
            plt.savefig(png_file, format="png", dpi=400)
        if not pdf_file and not png_file:
            plt.show()

        return self.accessions, lats, lons

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
        ax.legend(loc=2, numpoints=2, handlelen=0.05, markerscale=1, prop=fontProp, pad=0.018)
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
        ax.legend(loc=2, numpoints=2, handlelen=0.05, markerscale=1, prop=fontProp, pad=0.018)
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


