import holoviews as hv
import numpy as np
from holoviews import opts
hv.extension('bokeh')


def number_concentration_per_roi(d,
                                 section,
                                 version,
                                 iteration,
                                 rois_version,
                                 rois_per_class,
                                 ordering,
                                 signs=None,
                                 means=None,
                                 max_y=0):
    target = 'Total EVs'
    dt = []
    aggregateds = {}
    overlays = []
    labels = dict()
    for i in range(len(d['roi_areas'])):
        for j in range(len(d['total_evs_per_roi'][i])):
            # fetch the ROIs per class for this section
            for k, v in rois_per_class.items():
                if i + 1 in v:
                    # class, roi, count
                    dt.append((k, i, d['total_evs_per_roi'][i][j]))
                    # labels test
                    if signs is not None and i < len(signs):
                        if signs[i] < 0.05:
                            if k not in labels:
                                labels[k] = []
                            labels[k].append((i, means[i] * 1.1, '*'))

    dset = hv.Dataset(dt, ['Class', 'ROI', target], [target])
    for class_name, class_colour in ordering:
        aggregateds[class_name] = dset[class_name].aggregate(['ROI'], np.mean,
                                                             np.std)
        # if len(aggregateds[class_name][target]) > 0:
        # m = max(aggregateds[class_name][target])
        # if m > max_y:
        #    max_y = m
        if class_name in labels:
            overlays.append(
                hv.Bars(aggregateds[class_name]).opts(
                    color=class_colour, fill_alpha=0.5, title=class_name) *
                hv.ErrorBars(aggregateds[class_name]) *
                hv.Labels(labels[class_name]))
        else:
            overlays.append(
                hv.Bars(aggregateds[class_name]).opts(
                    color=class_colour, fill_alpha=0.5, title=class_name) *
                hv.ErrorBars(aggregateds[class_name]))

    # export this as an image
    bars = hv.Layout(
        overlays,
        group=f'Number concentration per ROI, {section}. Mean + SD',
        label='').opts(opts.Bars(ylim=(0, max_y * 1.1))).opts(toolbar=None)

    target_file = f"./resources/analysis/output/plot_{version}_{iteration}_rois_v{rois_version}_numconc_per_roi.png"
    print('number_concentration_per_roi saved in:')
    print(target_file)
    hv.save(bars, target_file, fmt='png', size=100)


def freq_bars(means1,
              sd1,
              means2,
              sd2,
              probs,
              ver1,
              iter1,
              ordering,
              rois_per_class,
              rois_ver1,
              ver2,
              iter2,
              rois_ver2,
              columns=6):
    plots = []
    maxy = max(np.max(means1), np.max(means2))
    maxsd = max(np.max(sd1), np.max(sd2))
    print('maxy', maxy, 'maxsd', maxsd)

    for cname, colour in ordering:
        # collect means, sds and labels for this class
        mns1, mns2 = [], []
        err1, err2 = [], []
        lbls1, lbls2 = [], []
        for i in range(len(means1)):
            if i + 1 in rois_per_class[cname]:
                mns1.append((f'{i+1}', means1[i]))
                mns2.append((f'{i+1}', means2[i]))
                err1.append((f'{i+1}', means1[i], sd1[i]))
                err2.append((f'{i+1}', means2[i], sd2[i]))
                if probs[i] < 0.05:
                    lbls1.append((f'{i+1}', means1[i] + sd1[i], '*'))
                    lbls2.append((f'{i+1}', means2[i] + sd2[i], '*'))

        if len(mns1):
            h1 = hv.Bars(mns1).options(
                color=colour,
                alpha=0.7,
                ylim=(0.1, (maxy + maxsd) * 1.01),
                title=f"Number concentration {cname}. Mean + SD. Pre-Ov",
                xlabel='EV radius in um',
                ylabel='Frequency',
                bgcolor="#E8DDCB")
            e1 = hv.ErrorBars(err1)

            h2 = hv.Bars(mns2).options(
                color=colour,
                alpha=0.7,
                ylim=(0.1, (maxy + maxsd) * 1.01),
                title=f"Number concentration {cname}. Mean + SD. Post-Ov",
                xlabel='EV radius in um',
                ylabel='Frequency',
                bgcolor="#E8DDCB")
            e2 = hv.ErrorBars(err2)

            if len(lbls1) > 0:
                plots.append(
                    (h1 * e1 * hv.Labels(lbls1) + h2 * e2 * hv.Labels(lbls2)))
            else:
                plots.append((h1 * e1 + h2 * e2))

    target_file = f"./resources/analysis/output/plot_[{ver1}_{iter1}_rois_v{rois_ver1}-{ver2}_{iter2}_rois_v{rois_ver2}]_numconc_per_roi.png"
    print('Size distribution per ROI saved in:')
    print(target_file)
    hv.save(hv.Layout(plots).cols(columns).opts(toolbar=None),
            target_file,
            fmt='png',
            size=200)


def size_distribution_per_roi_np(d,
                                 section,
                                 version,
                                 iteration,
                                 rois_per_class,
                                 replicates,
                                 ordering,
                                 bins=10,
                                 columns=2):
    dt = []
    # create storage for the counts per roi and repeat
    for roi in range(len(d['roi_areas'])):
        li = []
        for rep in range(replicates):
            li.append(list())
        dt.append(li)

    mi = 10
    ma = 0
    # find maximum and minimum values to use for the plotting range
    # row contains (roi, replicate, size)
    for row in d['evs_in_roi_replicate_radius_age']:
        if row[2] > ma:
            ma = row[2]
        elif row[2] < mi:
            mi = row[2]
        dt[row[0]][row[1]].append(row[2])

    plots = []
    maxy = 0

    # compute the values for the histograms and the error bars
    for roi in range(len(d['roi_areas'])):
        fs = []  # storage for the frequencies per plot
        for rep in range(replicates):
            freqs, edges = np.histogram(dt[roi][rep],
                                        bins=bins,
                                        range=(mi, ma))
            fs.append(freqs)
        # compute mean, sd per bin
        values = np.zeros([bins, replicates])
        for bi in range(bins):
            for rep in range(replicates):
                v = fs[rep][bi]
                maxy = v if v > maxy else maxy
                values[bi][rep] = fs[rep][bi]
        frequencies = values.mean(axis=1)
        devs = values.std(axis=1)

        for cname, color in ordering:
            if roi + 1 in rois_per_class[cname]:
                col = color
                break
        h = hv.Histogram(
            (edges, frequencies)).opts(color=col,
                                       fill_alpha=0.5,
                                       title=f"ROI #{roi + 1} size dist.",
                                       xlabel='EV radius in um',
                                       ylabel='Frequency',
                                       bgcolor="#E8DDCB",
                                       ylim=(0, maxy),
                                       padding=(0.05, 0))
        e = hv.ErrorBars((edges, frequencies, devs))
        plots.append((h * e))
    print('maxy:', maxy)
    target_file = f"./resources/analysis/output/plot_{version}_{iteration}_hist_sizedist_rois.png"
    print('Size distribution per ROI saved in:')
    print(target_file)
    hv.save(hv.Layout(plots).cols(columns).opts(toolbar=None),
            target_file,
            fmt='png',
            size=200)


def size_bars_per_roi(means1,
                      sd1,
                      means2,
                      sd2,
                      probs,
                      bin_edges,
                      ver1,
                      iter1,
                      ordering,
                      rois_per_class,
                      rois_ver1,
                      ver2,
                      iter2,
                      rois_ver2,
                      columns=6):
    plots = []
    bin_edges[0] = 0.04
    bin_edges[-1] = 0.60001
    maxy = max(np.max(means1), np.max(means2))
    maxsd = max(np.max(sd1), np.max(sd2))
    print('maxy', maxy, 'maxsd', maxsd)
    print('bin_edges', bin_edges)

    for cname, colour in ordering:
        for i in range(len(means1)):
            if i + 1 in rois_per_class[cname]:
                m1 = means1[i]
                m2 = means2[i]
                s1 = sd1[i]
                s2 = sd2[i]
                prbs = probs[i]
                labels1, labels2 = [], []
                print(m1, 'and', m2, prbs)

                h1 = hv.Histogram((bin_edges, m1)).options(
                    color=colour,
                    alpha=0.7,
                    ylim=(0.1, (maxy + maxsd) * 1.01),
                    title=f"ROI #{i + 1} size distribution. Mean + SD. Pre-Ov",
                    xlabel='EV radius in um',
                    ylabel='Frequency',
                    bgcolor="#E8DDCB")
                e1 = hv.ErrorBars((bin_edges, m1, s1))

                h2 = hv.Histogram((bin_edges, m2)).options(
                    color=colour,
                    alpha=0.7,
                    ylim=(0.1, (maxy + maxsd) * 1.01),
                    title=f"ROI #{i + 1} size distribution. Mean + SD. Post-Ov",
                    xlabel='EV radius in um',
                    ylabel='Frequency',
                    bgcolor="#C8DDCB")
                e2 = hv.ErrorBars((bin_edges, m2, s2))
                for j in range(len(m1)):
                    if prbs[j] < 0.05:
                        mp = (bin_edges[j] + bin_edges[j + 1]) / 2
                        labels1.append((mp, m1[j] + s1[j], '*'))
                        labels2.append((mp, m2[j] + s2[j], '*'))
                if len(labels1) > 0:
                    plots.append((h1 * e1 * hv.Labels(labels1) +
                                  h2 * e2 * hv.Labels(labels2)))
                else:
                    plots.append((h1 * e1 + h2 * e2))

    target_file = f"./resources/analysis/output/plot_[{ver1}_{iter1}_rois_v{rois_ver1}-{ver2}_{iter2}_rois_v{rois_ver2}]_hist_sizedist_rois.png"
    print('Size distribution per ROI saved in:')
    print(target_file)
    hv.save(hv.Layout(plots).cols(columns).opts(toolbar=None),
            target_file,
            fmt='png',
            size=200)
