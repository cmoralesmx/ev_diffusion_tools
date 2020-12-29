import bz2
import copy
import pickle
import numpy as np
import scipy.stats as sp_stats
from itertools import combinations
from ev_model import plotting  # , plotting_mpl

sections = ['isthmus', 'ampulla']

classes2 = {
    'ordering': [('narrow_end', '#30a2da'), ('wide_end', '#fc4f30'),
                 ('narrow_lumen', '#e5ae38'), ('wide_lumen', 'green')]
}
classes2['isthmus'] = {
    'narrow_end': [1, 2, 3, 4, 5, 6],
    'wide_end': [7, 8, 9],
    'narrow_lumen': [10, 11, 12, 13, 14, 15],
    'wide_lumen': [16, 17, 18, 19, 20, 21]
}
classes2['ampulla'] = {
    'narrow_end': [1, 2, 3, 4, 5, 6],
    'wide_end': [7, 8, 9, 10, 11, 12],
    'narrow_lumen': [13, 14, 15, 16, 17, 18],
    'wide_lumen': [19, 20, 21, 22, 23, 24]
}

# deffinitions for v3 of the analysis
# color palete from https://learnui.design/tools/data-color-picker.html#palette
classes3 = {
    'ordering': [('narrow_end', '#003f5c'), ('narrow_lumen', '#444e86'),
                 ('lumen_centre', '#955196'), ('wide_end', '#ffa600'),
                 ('near_epithelium', '#dd5182'), ('far_epithelium', '#ff6e54')]
}
classes3['isthmus'] = {
    'narrow_end': [1, 2, 4, 5, 6, 7],
    'wide_end': [],
    'narrow_lumen': [9, 10, 11, 12, 13, 14],
    'lumen_centre': [17, 18],
    'near_epithelium': [1, 2, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14],
    'far_epithelium': [3, 8, 15, 16, 19, 20, 21]
}
classes3['ampulla'] = {
    'narrow_end': [1, 2, 3, 4, 5, 6],
    'wide_end': [7, 8, 9, 10, 11, 12],
    'narrow_lumen': [13, 14, 15, 16, 17, 18],
    'lumen_centre': [19, 20, 21, 22, 23, 24],
    'near_epithelium': [32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43],
    'far_epithelium': [25, 26, 27, 28, 29, 30, 31]
}


def log_normality_test(data, log):
    """
    The result from D'Agostino and Pearson's test is reported
    """
    for p in data:
        log.write(
            f"\tnormality test {f'p={p:.3f}' if p > 1E-3 else f'p={p:.3E}'} ")
        if p < 0.05:
            log.write('(significantly non-normal)\n')
        else:
            log.write('(norm. dist.)\n')


def validate_parametric(min_n, dpr, alpha, log):
    parametric = False
    if min_n > 99:
        log.write(
            f'min_n ({min_n}) > 100, Normality not essential, relaxing this '
            f'assumption. Yet, still checking \n')
        # same variances?
        _, levene_p = sp_stats.levene(*dpr)
        parametric = True
        if levene_p > alpha:
            log.write('\tFail to reject Homoscedasticity\n')
        else:
            log.write(
                '\tSignificant differences in variance irrelevant (large N)\n')

    elif min_n > 25:
        log.write(f'25 < min_n ({min_n}) < 10, Normality test is required\n')
        tn = [sp_stats.normaltest(d)[1] for d in dpr]
        if min(tn) < alpha:
            log.write(
                '\n\tAt least one distribution is SIGNIFICANTLY non-normal:\n')
            log_normality_test(tn, log)
        else:
            log.write(
                '\n\tDistributions not significantly different from normal\n')
            # same variances?
            _, levene_p = sp_stats.levene(*dpr)
            if levene_p > alpha:
                parametric = True
                log.write('\tFail to reject Homoscedasticity\n')

            else:
                log.write('\tSIGNIFICANTLY DIFFERENT VARIANCES.\n')
    else:
        log.write(f'3 < min_n ({min_n}) < 25, \n')
    return parametric


def validate_double_parametric(min_n1, dpr1, min_n2, dpr2, alpha, log):
    parametric = False
    if min_n1 > 99 and min_n2 > 99:
        log.write(
            f'min_n1 ({min_n1}) > 100 and min_n2 ({min_n2}) > 99, Normality '
            f'not essential, relaxing this assumption. Yet, still checking\n')
        # same variances?
        _, levene_p = sp_stats.levene(*dpr1, *dpr2)
        parametric = True
        if levene_p > alpha:
            log.write('Fail to reject Homoscedasticity ')
        else:
            log.write('SIGNIFICANTLY DIFFERENT VARIANCES. ')

    elif min_n1 > 25 and min_n2 > 25:
        log.write(f'25 < min_n1 ({min_n1}), min_n2 ({min_n2}) < 100,\n')
        tn = [sp_stats.normaltest(d)[1] for d in dpr1 + dpr2]
        if min(tn) < alpha:
            log.write(
                '\n\tAt least one distribution is SIGNIFICANTLY non-normal. ')
        else:
            log.write(
                '\n\tDistributions not significantly different from normal. ')
            # same variances?
            _, levene_p = sp_stats.levene(*dpr1, *dpr2)
            if levene_p > alpha:
                parametric = True
                log.write('Fail to reject Homoscedasticity ')
            else:
                log.write('SIGNIFICANTLY DIFFERENT VARIANCES. ')
    else:
        log.write(f'3 < min_n ({min_n1}, {min_n2}) < 25, ')
    return parametric


def identify_testable_rois(dpr, lengths, log):
    min_n = min([len(d) for d in dpr])
    # samples with less than 3 data points cannot be tested
    if min_n < 3:
        log.write(
            f'min_n ({min_n}) < 3, Attempting testing a subset of the data. \n'
        )
        dpr = [d for d in dpr if len(d) > 3]
        if len(dpr) < 2:
            log.write('There are not sufficient data points for testing\n')
            return None, None
        log.write(
            f'Elements skipped: {[pair for pair in lengths if pair[1] < 3]} \n'
        )
        min_n = min([len(d) for d in dpr])
    return min_n, dpr


def frequencies_per_roi(data, bin_edges, rois, sec, n_rep, log):
    freqs_per_roi, binned_freq_per_roi = [], []
    binned_freq_per_roi_mean, binned_freq_per_roi_sd = [], []
    for i in range(len(rois)):
        freqs, binned_freqs, = [], []
        for rep in range(n_rep):
            d = data.query(
                f"section=='{sec}' & roi=={rois[i]} & replicate=={rep}"
            )['radius_um']
            # frequency per ROI in this repeat
            freqs.append(len(d))
            # binned frequency per ROI
            bf, _ = np.histogram(d, bins=bin_edges)
            binned_freqs.append(bf)

        log.write(f'ROI #{rois[i]}, Total EV per repeat: {freqs}\n')
        log.write(
            f'Total EV p/size-range p/repeat(binned_freqs) {binned_freqs}\n')

        freqs_per_roi.append(freqs)
        binned_freq_per_roi.append(binned_freqs)

        npa = np.array(binned_freqs)
        log.write(f'Mean EV per size-range: {npa.mean(axis=0)}\n')
        binned_freq_per_roi_mean.append(npa.mean(axis=0))
        binned_freq_per_roi_sd.append(npa.std(axis=0))

    log.write('Summary statistics for cross-section:\n')
    freqs_per_roi_means = np.array(freqs_per_roi).mean(axis=1)
    freqs_per_roi_sd = np.array(freqs_per_roi).std(axis=1)
    freqs_per_roi_var = np.array(freqs_per_roi).var(axis=1)
    log.write(f'Mean EV { freqs_per_roi_means }\n')
    log.write(f'Std Dev { freqs_per_roi_sd }\n')
    log.write(f'Varianc { freqs_per_roi_var }\n')
    return [
        freqs_per_roi, freqs_per_roi_means, freqs_per_roi_sd,
        freqs_per_roi_var, binned_freq_per_roi, binned_freq_per_roi_mean,
        binned_freq_per_roi_sd
    ]


def test_groups(dpr1, dpr2, alpha, parametric, log):
    if parametric:
        res = sp_stats.f_oneway(*dpr1, *dpr2)
        print(f'\n\tParametric test (one-way ANOVA). {res[1]:.4f}', end=' ')
        log.write(f'\n\tParametric test (one-way ANOVA). {res[1]:.4f} ')
    else:
        res = sp_stats.kruskal(*dpr1, *dpr2)
        print(f'\n\tNon-parametric test (Kruskal-Wallis). {res[1]:.4f}',
              end=' ')
        log.write(f'\n\tNon-parametric test (Kruskal-Wallis). {res[1]:.4f} ')
    print('NON sig. diffs'
          if res[1] > alpha else 'SIGNIFICANT DIFFERENCES. Reject H0')
    log.write('NON sig. diffs'
              if res[1] > alpha else 'SIGNIFICANT DIFFERENCES. Reject H0')
    dprs = dpr1 + dpr2
    pooled_dev = np.sqrt(
        np.sum([(r.count() - 1) * (r.mean()**2)
                for r in dprs]) / (np.sum([r.count()
                                           for r in dprs]) - len(dprs)))
    print(f'\n    Pooled deviation {pooled_dev}\n\n')
    log.write(f'\n    Pooled deviation {pooled_dev}\n\n')


def test_group(dpr, alpha, parametric, log):
    if parametric:
        res = sp_stats.f_oneway(*dpr)
        print(f'\n\tParametric test (one-way ANOVA). {res[1]:.4f}', end=' ')
        log.write(f'\n\tParametric test (one-way ANOVA). {res[1]:.4f} ')
    else:
        res = sp_stats.kruskal(*dpr)
        print(f'\n\tNon-parametric test (Kruskal-Wallis). {res[1]:.4f}',
              end=' ')
        log.write(f'\n\tNon-parametric test (Kruskal-Wallis). {res[1]:.4f} ')

    if res[1] > alpha:
        print('NON sig. diffs')
        log.write('NON sig. diffs')

        pooled_dev = np.sqrt(
            np.sum([(r.count() - 1) * (r.mean()**2)
                    for r in dpr]) / (np.sum([r.count()
                                              for r in dpr]) - len(dpr)))
        print(f'\n    Pooled deviation {pooled_dev}\n\n')
        log.write(f'\n    Pooled deviation {pooled_dev}\n\n')
    else:
        print('SIGNIFICANT DIFFERENCES. Reject H0')
        log.write('SIGNIFICANT DIFFERENCES. Reject H0\n')
        # extra test should be carried to identify pairwise differences


def export_csv(sections, versions, iter1, rois1, freqs_per_roi1,
               freqs_per_roi2, binned_freqs_per_roi1, binned_freqs_per_roi2,
               rois_per_class, apop, differentiation):
    s1, v1 = sections[0], versions[0]
    s2, v2 = sections[1], versions[1]

    class_per_roi = []
    drois = ['near_epithelium', 'far_epithelium']
    # produce a list of tuples [(roi_id, class)]
    for c, rois in rois_per_class.items():
        for r in rois:
            class_per_roi.append((r - 1, c))

    f_rois_p = f'./resources/analysis/output/freqs_rois_{s1}_{v1}-'\
               f'{s2}_{v2}-{iter1}_rois_v{rois1}.csv'

    f_size_p = f'./resources/analysis/output/freqs_size_{s1}_{v1}-{s2}'\
               f'_{v2}-{iter1}_rois_v{rois1}.csv'
    if rois1 == '4.0':
        f_distance_p = f'./resources/analysis/output/freqs_distance_{s1}_{v1}'\
                       f'-{s2}_{v2}-{iter1}_rois_v{rois1}.csv'
    else:
        f_distance_p = "/dev/null"

    print('CSV freqs per ROI, saved as:')
    print(f_rois_p)
    print(f_size_p)
    if rois1 == '4.0':
        print(f_distance_p)

    if iter1 == 14400000:
        t = '4'
    elif iter1 == 28800000:
        t = '8'
    else:
        t = '10'

    sec = s1[:-1] if s1[-1:] == '1' else s1

    # needs info about the roi IDs
    with open(f_rois_p, 'w') as f_rois:
        f_rois.write(
            'section,stage,roi,roi_type,time,freq,apoptosis,differentiation\n')
        n_rois = len(freqs_per_roi1)
        for (freqs_per_roi, stage) in [(freqs_per_roi1, 'pre'),
                                       (freqs_per_roi2, 'post')]:
            for roi_id in range(n_rois):
                # fin the type of ROI
                n_repeats = len(freqs_per_roi[roi_id])
                for repeat_id in range(n_repeats):
                    freq = freqs_per_roi[roi_id][repeat_id]
                    f_rois.write(
                        f'{sec},{stage},{class_per_roi[roi_id][0]},'
                        f'{class_per_roi[roi_id][1]},{t},{freq},{apop},'
                        f'{differentiation}\n')

    with open(f_size_p, 'w') as f_size, open(f_distance_p, 'w') as f_distance:
        f_size.write('section,stage,roi,roi_type,time,freq,apoptosis,'
                     'differentiation,type\n')
        if rois1 == '4.0':
            f_distance.write('section,stage,roi,roi_type,time,freq,apoptosis,'
                             'differentiation,type\n')

        n_rois = len(freqs_per_roi1)
        for (binned_freqs_per_roi, stage) in [(binned_freqs_per_roi1, 'pre'),
                                              (binned_freqs_per_roi2, 'post')]:
            for roi_id in range(n_rois):
                n_repeats = len(binned_freqs_per_roi[roi_id])
                for repeat_id in range(n_repeats):
                    mvs = 0
                    for i_size in range(
                            len(binned_freqs_per_roi[roi_id][repeat_id])):
                        if i_size < 1:
                            exo = binned_freqs_per_roi[roi_id][repeat_id][
                                i_size]
                            msg = f'{sec},{stage},' \
                                  f'{class_per_roi[roi_id][0]},' \
                                  f'{class_per_roi[roi_id][1]},' \
                                  f'{t},{exo},{apop},' \
                                  f'{differentiation},exosomes\n'
                            if class_per_roi[roi_id][1] in drois:
                                f_distance.write(msg)
                            else:
                                f_size.write(msg)
                        else:
                            mvs += binned_freqs_per_roi[roi_id][repeat_id][
                                i_size]
                    # write the values for microvesicles
                    msg = f'{sec},{stage},{class_per_roi[roi_id][0]},' \
                          f'{class_per_roi[roi_id][1]},{t},{mvs},{apop},' \
                          f'{differentiation},microvesicles\n'
                    if class_per_roi[roi_id][1] in drois:
                        f_distance.write(msg)
                    else:
                        f_size.write(msg)


def differences_across_classes(log,
                               data,
                               section,
                               replicates,
                               pooled_repeats=False,
                               alpha=0.05):
    """
    Compare the ROIs in one class with those in the other classes.
    Example:
        From the ROIS in the isthmus, the ROIs classed as 'narrow-end'
        can be compared agains the ROIs classed as 'wide_lumen',
        'wide-end', and 'narrow-lumen'
    """
    li = f"{'=='*20}\nComparing ROIs across classes in {section} "\
        f"{'POOLING REPEATS' if pooled_repeats else ''}\n"
    print(li)
    print('H0: the samples were drawn from a pop. with the same distribution')
    log.write(li)
    log.write('H0: the samples/group of samples were drawn from a population'
              'with the same distribution\n')

    pairs_of_classes = [
        p for p in combinations([t[0] for t in classes['ordering']], 2)
    ]

    for roi_class_1, roi_class_2 in pairs_of_classes:
        rois_class1 = [roi for roi in classes[section][roi_class_1]]
        rois_class2 = [roi for roi in classes[section][roi_class_2]]
        print(f"{'=='*10}\n")
        print('Rois in', roi_class_1, rois_class1)
        print('Rois in', roi_class_2, rois_class2)
        log.write(f"{'=='*10}\n")
        log.write(f'Rois in {roi_class_1} {rois_class1}\n')
        log.write(f'Rois in {roi_class_2} {rois_class2}\n')

        for rep in range(1 if pooled_repeats else replicates):
            if not pooled_repeats:
                print(f'------------ repeat {rep+1} ------------')
                log.write(f'------------ repeat {rep+1} ------------\n')

            process_repeat(data, rep, section, section, replicates, replicates,
                           rois_class1, rois_class2, log, alpha)


# comparisons for two cross sections
def differences_two_sections(log,
                             data,
                             sections,
                             versions,
                             replicates,
                             bin_edges,
                             classes,
                             paired_classes=True,
                             alpha=0.05):
    """
    Compares the ROIs of a given class in the first section against
    ROIs of non-matching classes in the second section.
    This comparison is based on comparing the mean values after binning the
    EVs per size
    """
    mc = ('in MATCHING classes in'
          if paired_classes else 'in NON-matching classes in')
    li = f"{'=='*20}\nComparing ROIs {mc} {sections}\n"
    print(li)
    print('H0: both samples were drawn from a pop. with the same distribution')
    log.write(li)
    log.write('H0: both samples/group of samples were drawn from a population '
              'with the same distribution\n')

    class_reference = ['', '']
    if type(sections) != list:
        raise RuntimeError('sections is not a list!', type(sections))
    else:
        # remove the added 1 and 2 numbers if needed
        log.write(f'Renaming sections {sections}\n')
        class_reference[
            0] = sections[0][:-1] if sections[0][-1] == '1' else sections[0]
        class_reference[
            1] = sections[1][:-1] if sections[1][-1] == '2' else sections[1]
        log.write(f'New names {class_reference}\n')
    if paired_classes:
        pairs_of_classes = [(t[0], t[0]) for t in classes['ordering']]
    else:
        pairs_of_classes = [
            p for p in combinations([t[0] for t in classes['ordering']], 2)
        ]

    log.write(f'pairs_of_classes: {pairs_of_classes}\n')
    n_rois1, n_rois2 = 0, 0
    # create a list of ROIs
    rois_class1, rois_class2 = [], []
    for class_1, class_2 in pairs_of_classes:
        rois_class1 += [
            i_roi - 1 for i_roi in classes[class_reference[0]][class_1]
        ]
        rois_class2 += [
            i_roi - 1 for i_roi in classes[class_reference[1]][class_2]
        ]
        n_rois1 += len(rois_class1)
        n_rois2 += len(rois_class2)
    log.write(f'rois_class1: {rois_class1}\n')
    log.write(f'rois_class2: {rois_class2}\n')

    log.write('=== Frequencies for section 1\n')
    [
        freqs_per_roi1, freqs_per_roi_means1, freqs_per_roi_sd1,
        freqs_per_roi_var1, binned_freq_per_roi1, binned_freq_per_roi_mean1,
        binned_freq_per_roi_sd1
    ] = frequencies_per_roi(data, bin_edges, rois_class1, sections[0],
                            replicates[0], log)
    log.write('=== Frequencies for section 2\n')
    [
        freqs_per_roi2, freqs_per_roi_means2, freqs_per_roi_sd2,
        freqs_per_roi_var2, binned_freq_per_roi2, binned_freq_per_roi_mean2,
        binned_freq_per_roi_sd2
    ] = frequencies_per_roi(data, bin_edges, rois_class2, sections[1],
                            replicates[1], log)

    max_y = (max(np.max(freqs_per_roi_means1), np.max(freqs_per_roi_means2)) +
             max(np.max(freqs_per_roi_sd1), np.max(freqs_per_roi_sd2)))

    return [
        freqs_per_roi1, freqs_per_roi2, freqs_per_roi_means1,
        freqs_per_roi_sd1, freqs_per_roi_means2, freqs_per_roi_sd2,
        binned_freq_per_roi1, binned_freq_per_roi2, binned_freq_per_roi_mean1,
        binned_freq_per_roi_sd1, binned_freq_per_roi_mean2,
        binned_freq_per_roi_sd2, max_y
    ]


def header_double(log, section1, version1, path1, iter1, rois1, rep1, section2,
                  version2, path2, iter2, rois2, rep2, paired_classes,
                  number_concentration, size_distribution, bins, columns):
    print('Comparing->        [1]   |   [2]')
    print(f'   section:      {section1} | {section2}')
    print(f'   version:  {version1} | {version2}')
    print(
        f' iteration:        {int(iter1)/1e6:>.1f}M | {int(iter2)/1e6:>.1f}M')
    print(f'      rois:           {rois1} | {rois2}')
    print(f'replicates:            {rep1} | {rep2}')
    print(f'path [1]: {path1}')
    print(f'path [2]: {path2}')
    print('------------------------------------')
    print('| paired classes?', paired_classes)
    print('| plot number concentration?', number_concentration)
    print('| plot size distribution?', size_distribution)
    print('| number of bins:', bins, 'number of columns:', columns)
    print('------------------------------------')
    log.write('Comparing->        [1]   |   [2]\n')
    log.write(f'   section:      {section1} | {section2}\n')
    log.write(f'   version:  {version1} | {version2}\n')
    log.write(
        f' iteration:        {int(iter1)/1e6:>.1f}M | {int(iter2)/1e6:>.1f}M\n'
    )
    log.write(f'      rois:           {rois1} | {rois2}\n')
    log.write(f'replicates:            {rep1} | {rep2}\n')
    log.write(f'path [1]: {path1}\n')
    log.write(f'path [2]: {path2}\n')
    log.write('------------------------------------\n')
    log.write(f'| paired classes? {paired_classes}\n')
    log.write(f'| plot number concentration? {number_concentration}\n')
    log.write(f'| plot size distribution? {size_distribution}\n')
    log.write(f'| number of bins: {bins}, number of columns: {columns}\n')
    log.write('------------------------------------\n')


def compare_double(log,
                   section1,
                   version1,
                   path1,
                   iter1,
                   rois1,
                   rep1,
                   section2,
                   version2,
                   path2,
                   iter2,
                   rois2,
                   rep2,
                   number_concentration,
                   size_distribution,
                   bins,
                   columns,
                   apoptosis,
                   controlled,
                   paired_classes=True):

    header_double(log, section1, version1, path1, iter1, rois1, rep1, section2,
                  version2, path2, iter2, rois2, rep2, paired_classes,
                  number_concentration, size_distribution, bins, columns)
    stats = []
    data, data_per_replicate = None, None
    source1 = f'./resources/analysis/output/stats_{version1}_{iter1}'\
              f'_rois_v{rois1}.pickle.bz2'
    source2 = f'./resources/analysis/output/stats_{version2}_{iter2}'\
              f'_rois_v{rois2}.pickle.bz2'
    for source in [source1, source2]:
        with bz2.BZ2File(source, 'rb') as stats_source:
            s = pickle.load(stats_source)

            df = s['data_frame']
            df_per_replicate = s['evs_per_replicate_df']

            if data is None:
                if section1 == section2:
                    print('adding 1 to the section name')
                    df.index = df.index.set_levels(df.index.levels[0].str.cat(
                        ['1']),
                                                   level=0)
                    df_per_replicate.index = df_per_replicate.index.set_levels(
                        df_per_replicate.index.levels[0].str.cat(['1']),
                        level=0)
                data = df
                data_per_replicate = df_per_replicate
            else:
                if section1 == section2:
                    print('adding 2 to the section name')
                    df.index = df.index.set_levels(df.index.levels[0].str.cat(
                        ['2']),
                                                   level=0)
                    df_per_replicate.index = df_per_replicate.index.set_levels(
                        df_per_replicate.index.levels[0].str.cat(['2']),
                        level=0)
                data = data.append(s['data_frame'])
                data_per_replicate = data_per_replicate.append(
                    s['evs_per_replicate_df'])
            stats.append(s)
    if section1 == section2:
        s1 = section1 + '1'
        s2 = section2 + '2'
    else:
        s1 = section1
        s2 = section2

    # produce mean and sd of exosomes and micro vesicles per replicate
    exs1, exs2, mvs1, mvs2, t1, t2 = [], [], [], [], [], []
    e17 = '=' * 17
    m15 = '-' * 15
    m16 = '-' * 16
    m36 = '-' * 36
    print('\n', e17, 'Frequency of Exosomes and MVs per repeat', e17)
    print(' ', m15, [1], m16, ' ', m15, [2], m16)
    print('  |   Exos    |    MVs    ||    Sum    |'
          ' |   Exos    |    MVs    ||    Sum    |')
    print(f'   {m36}| |{m36}|')
    for r in range(rep1):
        e1 = data_per_replicate.query(
            f"section=='{s1}' & replicate=={r} & radius_um < 0.1"
        )['radius_um'].count()
        m1 = data_per_replicate.query(
            f"section=='{s1}' & replicate=={r} & radius_um >= 0.1"
        )['radius_um'].count()
        e2 = data_per_replicate.query(
            f"section=='{s2}' & replicate=={r} & radius_um < 0.1"
        )['radius_um'].count()
        m2 = data_per_replicate.query(
            f"section=='{s2}' & replicate=={r} & radius_um >= 0.1"
        )['radius_um'].count()

        print(f'  | {e1:>6}    | {m1:>6}    || {e1 + m1:>6}    ', end='| |')
        print(f' {e2:>6}    | {m2:>6}    || {e2 + m2:>6}    |')
        exs1.append(e1)
        mvs1.append(m1)
        t1.append(e1 + m1)
        exs2.append(e2)
        mvs2.append(m2)
        t2.append(e2 + m2)
    print(f'{m36}---| |{m36}|')
    print(
        f'μ | {np.array(exs1).mean():>9.2f} | {np.array(mvs1).mean():>9.2f} || {np.array(t1).mean():>9.2f} |',
        end='')
    print(
        f' | {np.array(exs2).mean():>9.2f} | {np.array(mvs2).mean():>9.2f} || {np.array(t2).mean():>9.2f} |'
    )
    print(
        f'σ | {np.array(exs1).std():>9.2f} | {np.array(mvs1).std():>9.2f} || {np.array(t1).std():>9.2f} |',
        end='')
    print(
        f' | {np.array(exs2).std():>9.2f} | {np.array(mvs2).std():>9.2f} || {np.array(t2).std():>9.2f} |'
    )
    print(
        '---------------------------------------   ------------------------------------'
    )

    print('\n\nFrom ROIS')
    c1 = data.query(f"section=='{s1}' & replicate==1")['radius_um'].count()
    c2 = data.query(f"section=='{s2}' & replicate==1")['radius_um'].count()
    print(f'EVs in ROIS in replicate 1 from [1]: {c1}')
    print(f'EVs in ROIS in replicate 1 from [2]: {c2}')
    log.write(f'EVs in ROIS in replicate 1 from [1]: {c1}\n')
    log.write(f'EVs in ROIS in replicate 1 from [2]: {c2}\n')
    if rois1 == '3.0':
        print('ROIS v3.0')
        classes = copy.deepcopy(classes2)
        bin_edges = [
            0, 0.090000001, 0.192000001, 0.2940000001, 0.39600001, 0.4980001,
            0.6
        ]
    else:
        print('ROIS v4.0')
        classes = copy.deepcopy(classes3)
        bin_edges = [
            0, 0.090000001, 0.192000001, 0.2940000001, 0.39600001, 0.4980001,
            0.6
        ]

    [
        freqs_per_roi1, freqs_per_roi2, means1, sds1, means2, sds2,
        binned_freq_per_roi1, binned_freq_per_roi2, binned_freq_per_roi_mean1,
        binned_freq_per_roi_sd1, binned_freq_per_roi_mean2,
        binned_freq_per_roi_sd2, max_y
    ] = differences_two_sections(log,
                                 data, [s1, s2], [version1, version2],
                                 [rep1, rep2],
                                 bin_edges,
                                 classes,
                                 paired_classes=paired_classes)

    # exporting data in table-like format for 3rd party analysis
    export_csv([s1, s2], [version1, version2], iter1, rois1, freqs_per_roi1,
               freqs_per_roi2, binned_freq_per_roi1, binned_freq_per_roi2,
               classes[section1], apoptosis, controlled)

    # compare frequencies per ROI
    p_rois = []
    for freqs1, freqs2 in zip(freqs_per_roi1, freqs_per_roi2):
        _, pvalue = sp_stats.f_oneway(freqs1, freqs2)
        log.write(f'Comparing {freqs1} and {freqs2} p={pvalue:.4f}\n')
        p_rois.append(pvalue)

    # compare binned frequencies per ROI
    p_rois_binned = []
    for rois1, rois2 in zip(binned_freq_per_roi1, binned_freq_per_roi2):
        pbins = []
        for b1, b2 in zip(zip(*rois1), zip(*rois2)):
            _, pvalue = sp_stats.f_oneway(b1, b2)
            log.write(f'Comparing {b1} and {b2} p={pvalue:.4f}\n')
            pbins.append(pvalue)
        p_rois_binned.append(pbins)

    # produce plots using the data in stats which comes with the dataframes
    # number concentration per ROI
    if number_concentration:
        plotting.freq_bars(means1, sds1, means2, sds2, p_rois, version1, iter1,
                           classes['ordering'], classes[section1], rois1,
                           version2, iter2, rois2)

    # size distribution per ROI
    if size_distribution:
        plotting.size_bars_per_roi(
            binned_freq_per_roi_mean1, binned_freq_per_roi_sd1,
            binned_freq_per_roi_mean2, binned_freq_per_roi_sd2, p_rois_binned,
            bin_edges, version1, iter1, classes['ordering'], classes[section1],
            rois1, version2, iter2, rois2)


if __name__ == '__main__':
    import argparse
    from datetime import datetime
    """
    comparisons possible:
    1 section - 1) within classes, 2) across classes
    2 sections - 1) within classes, 2) across classes
    """

    parser = argparse.ArgumentParser(add_help=True)
    s1 = parser.add_argument_group(
        'section1', 'Settings for the first section to compare')
    s1.add_argument('section',
                    choices=sections,
                    help='The 1st section to analyse.')
    s1.add_argument(
        'version',
        help="Input data version to read. Forms part of the file name to"
        "load: 'stats_SECTION_VERSION_ROIS.pickle.bz2 ")
    s1.add_argument('path',
                    metavar='base_path',
                    help='Directory where the target files are located '
                    '(1st section to compare)')
    s1.add_argument('rois',
                    help='Version of the ROIs used for producing the stats')
    s1.add_argument('apoptosis', help='Apoptosis time used during simulatiion')
    s1.add_argument(
        '--iteration1',
        nargs='?',
        type=int,
        help='Iteration number to read - 1st section. Default 14.4k')
    s1.add_argument(
        '--replicates1',
        nargs='?',
        type=int,
        help='The number of replicates to read - 1st section. Default 3')
    s1.add_argument(
        '--plotNC',
        action='store_true',
        help='Enables plotting the number concentration p/section and p/ROI')
    s1.add_argument(
        '--plotSD',
        action='store_true',
        help='Enables plotting the size distribution per section and per ROI')
    s1.add_argument('--columns',
                    nargs='?',
                    type=int,
                    help='Specify the number of columns of plots. Default: 3')
    s1.add_argument(
        '--bins',
        nargs='?',
        type=int,
        help='Specify the number of bins used for comparison. Default: 15')
    s1.add_argument('--minRadius', type=float, default=0.04)
    s1.add_argument('--maxRadius', type=float, default=0.5)
    s1.add_argument('--controlled',
                    action='store_true',
                    help='Specifies the distribution of the secretory cells is'
                    'the same in the repeats')

    s2 = parser.add_argument_group(
        'section2', """Settings for a second section to compare.
        Section2 is optional, however, it's required arguments are:
        --other_section, --other_roi_class, --other_path
        """)
    s2.add_argument('--section2',
                    choices=sections,
                    help='The 2nd section to analyse.')
    s2.add_argument('--version2',
                    help="Input data version to read."
                    "Forms part of the file name to load:i"
                    " 'stats_SECTION_VERSION.pickle.bz2 ")
    s2.add_argument('--path2',
                    help='Directory where the target files are '
                    'located (2nd section to compare)')
    s2.add_argument('--rois2',
                    help='Version of the ROIs used for producing the '
                    'stats (2nd section to compare)')
    s2.add_argument(
        '--iteration2',
        nargs='?',
        type=int,
        help='Iteration number to read - 2nd section. Default 14.4k')
    s2.add_argument(
        '--replicates2',
        nargs='?',
        type=int,
        help='The number of replicates to read - 2nd section. Default 3')
    s2.add_argument('--pairedClasses',
                    action='store_true',
                    help='If present, only ROIs from both sections in '
                    'matching classes are compared')
    args = parser.parse_args()

    iter1 = args.iteration1 if args.iteration1 else 14400000
    i1 = f'{iter1/1e6:.1f}'
    r1 = args.rois
    rep1 = args.replicates1 if args.replicates1 else 3
    columns = args.columns if args.columns else 3
    bins = args.bins if args.bins else 15

    tstamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    s, v = args.section, args.version

    iter2 = args.iteration2 if args.iteration2 else 14400000
    i2 = f'{iter2/1e6:.1f}'
    r2 = args.rois2
    rep2 = args.replicates2 if args.replicates2 else 3
    if not args.version2:
        raise RuntimeError('version2 not provided')

    pc = '-pairedClasses' if args.pairedClasses else '-nonPairedClasses'
    s2 = args.section2
    v2 = args.version2
    target_log = f"./resources/analysis/output/log_analysis_[{v}_{i1}"\
                 f"M_{r1}-{v2}_{i2}M_{r2}]{pc}_{tstamp}.txt"
    with open(target_log, 'w') as log:
        compare_double(log, args.section, args.version, args.path, iter1, r1,
                       rep1, args.section2, args.version2, args.path2, iter2,
                       r2, rep2, args.plotNC, args.plotSD, bins, columns,
                       args.apoptosis,
                       'controlled' if args.controlled else 'std')
    print('Log saved to:', target_log)
