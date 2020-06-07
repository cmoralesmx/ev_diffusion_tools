import bz2
import pickle
import numpy as np
import scipy.stats as sp_stats
from itertools import product, combinations
from ev_model import plotting #, plotting_mpl

sections = ['isthmus', 'ampulla']

classes2 = {'ordering':[('narrow_end','#30a2da'), ('wide_end','#fc4f30'),
            ('narrow_lumen','#e5ae38'),('wide_lumen','green')]}
classes2['isthmus'] = {'narrow_end':[1,2,3,4,5,6], 
                        'wide_end':[7,8, 9],
                        'narrow_lumen':[10, 11, 12, 13, 14, 15],
                        'wide_lumen':[16,17,18,19,20,21]}
classes2['ampulla'] = {'narrow_end':[1,2,3,4,5,6],
                        'wide_end':[7,8,9,10, 11, 12],
                        'narrow_lumen':[13,14,15,16,17,18], 
                        'wide_lumen':[19, 20, 21, 22, 23, 24]}

def log_normality_test(data, log):
    """
    The result from D'Agostino and Pearson's test is reported
    """
    for p in data:    
        log.write(f"\tnormality test {f'p={p:.3f}' if p > 1E-3 else f'p={p:.3E}'} ")
        if p < 0.05:
            log.write('(significantly non-normal)\n')
        else:
            log.write('(norm. dist.)\n')

def validate_parametric(min_n, dpr, alpha, log):
    parametric = False
    if min_n > 99:
        log.write(f'min_n ({min_n}) > 100, Normality not essential, relaxing this assumption. Yet, still checking \n')
        # same variances?
        _, levene_p = sp_stats.levene(*dpr)
        parametric = True
        if levene_p > alpha:
            log.write('\tFail to reject Homoscedasticity\n')
        else:
            log.write('\tSignificant differences in variance not relevant due to a large N\n')

    elif min_n > 25:
        log.write(f'25 < min_n ({min_n}) < 10, Normality test is required\n')
        tn = [sp_stats.normaltest(d)[1] for d in dpr]
        if min(tn) < alpha:
            log.write('\n\tAt least one distribution is SIGNIFICANTLY non-normal:\n')
            log_normality_test(tn, log)
        else:
            log.write('\n\tDistributions not significantly different from normal\n')
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
        log.write(f'min_n1 ({min_n1}) > 100 and min_n2 ({min_n2}) > 99, Normality not essential, relaxing this assumption. Yet, still checking\n')
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
            log.write('\n\tAt least one distribution is SIGNIFICANTLY non-normal. ')
        else:
            log.write('\n\tDistributions not significantly different from normal. ')
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
        log.write(f'min_n ({min_n}) < 3, Attempting testing a subset of the data. \n')
        dpr = [d for d in dpr if len(d)>3]
        if len(dpr) < 2:
            log.write('There are not sufficient data points for testing\n')
            return None, None
        log.write(f'Elements skipped: {[pair for pair in lengths if pair[1] < 3]} \n      ')
        min_n = min([len(d) for d in dpr])
    return min_n, dpr

# single cross-section comparisons
def differences_within_classes(log, data, section, replicates, pooled_repeats=False, alpha=0.05):
    """
    Compare the ROIs in a class with each other, i.e: All ROIs in the class
    'narrow-end' are compared in a single comparison.
    """
    l = f"{'=='*20}'\nComparing ROIs within classes in {section} {'POOLING REPEATS' if pooled_repeats else ''}\n"
    print(l)
    print('H0: the samples were drawn from a population with the same distribution')
    log.write(l)
    log.write('H0: both samples/group of samples were drawn from a population with the same distribution\n')

    for class_name, _ in classes2['ordering']:       
        rois_class1 = [roi for roi in classes2[section][class_name]]
        log.write(f"{'=='*10}\n")
        log.write(f'Rois in {class_name} {rois_class1}\n')
        print(f"{'=='*10}\n")
        print(f'ROIs: in {class_name} {rois_class1}\n')
    
        for rep in range(1 if pooled_repeats else replicates):
            if not pooled_repeats:
                print(f'-------- repeat: {rep + 1} --------')
                log.write(f'-------- repeat: {rep + 1} --------\n')
            
            r = replicates if pooled_repeats else rep
            dpr, lengths = select_data(data, rois_class1, section, r, pooled_repeats, log)

            min_n, dpr = identify_testable_rois(dpr, lengths, log)
            if min_n is None or dpr is None:
                continue
            
            # one way t-test assumes homoscedasticity, we must check the samples have the same variance
            parametric = validate_parametric(min_n, dpr, alpha, log)

            # do the test
            test_group(dpr, alpha, parametric, log)

def select_data(data, rois_class, sec, n_rep, pooled_repeats, log):
    dpr, lengths = [], []
    rp = f' in {[i for i in range(n_rep)]}' if pooled_repeats else f'=={n_rep}'
    for i_roi in rois_class:
        #rs = [i for i in range(n_rep)]
        d = data.query(f"section=='{sec}' & roi=={i_roi-1} & replicate{rp}")['radius_um']
        lengths.append((i_roi, len(d)))
        log.write(f'roi {i_roi}, elements: {len(d)}\n')
        dpr.append(d)
    return dpr, lengths

def test_groups(dpr1, dpr2, alpha, parametric, log):
    if parametric:
        res = sp_stats.f_oneway(*dpr1, *dpr2)
        print(f'\n\tParametric test (one-way ANOVA). {res[1]:.4f}', end=' ')
        log.write(f'\n\tParametric test (one-way ANOVA). {res[1]:.4f} ')
    else:
        res = sp_stats.kruskal(*dpr1, *dpr2)
        print(f'\n\tNon-parametric test (Kruskal-Wallis). {res[1]:.4f}', end=' ')
        log.write(f'\n\tNon-parametric test (Kruskal-Wallis). {res[1]:.4f} ')
    print('NON sig. diffs' if res[1] > alpha else 'SIGNIFICANT DIFFERENCES. Reject H0')
    log.write('NON sig. diffs' if res[1] > alpha else 'SIGNIFICANT DIFFERENCES. Reject H0')
    dprs = dpr1 + dpr2
    pooled_dev = np.sqrt(np.sum([(r.count()-1) * (r.mean()**2) for r in dprs]) / (np.sum([r.count() for r in dprs]) - len(dprs)))
    print(f'\n    Pooled deviation {pooled_dev}\n\n')
    log.write(f'\n    Pooled deviation {pooled_dev}\n\n')

def test_group(dpr, alpha, parametric, log):
    if parametric:
        res = sp_stats.f_oneway(*dpr)
        print(f'\n\tParametric test (one-way ANOVA). {res[1]:.4f}', end=' ')
        log.write(f'\n\tParametric test (one-way ANOVA). {res[1]:.4f} ')
    else:
        res = sp_stats.kruskal(*dpr)
        print(f'\n\tNon-parametric test (Kruskal-Wallis). {res[1]:.4f}', end=' ')
        log.write(f'\n\tNon-parametric test (Kruskal-Wallis). {res[1]:.4f} ')

    if res[1] > alpha:
        print('NON sig. diffs')
        log.write('NON sig. diffs')

        pooled_dev = np.sqrt(np.sum([(r.count()-1) * (r.mean()**2) for r in dpr]) / (np.sum([r.count() for r in dpr]) - len(dpr)))
        print(f'\n    Pooled deviation {pooled_dev}\n\n')
        log.write(f'\n    Pooled deviation {pooled_dev}\n\n')
    else:
        print('SIGNIFICANT DIFFERENCES. Reject H0')
        log.write('SIGNIFICANT DIFFERENCES. Reject H0\n')
        # extra test should be carried to identify pairwise differences

def process_repeat(data, rep, pooled_repeats, sec1, sec2, n_reps1, n_reps2, rois_class1, rois_class2, log, alpha):
    r1 = n_reps1 if pooled_repeats else rep
    r2 = n_reps2 if pooled_repeats else rep
    dpr1, lengths1 = select_data(data, rois_class1, sec1, r1, pooled_repeats, log)
    dpr2, lengths2 = select_data(data, rois_class2, sec2, r2, pooled_repeats, log)
    # all the test must be done using both lists of data
    min_n1, dpr1 = identify_testable_rois(dpr1, lengths1, log)
    min_n2, dpr2 = identify_testable_rois(dpr2, lengths2, log)
    
    if min_n1 is None or min_n2 is None:
        return None

    # one way t-test assumes homoscedasticity, we must check the samples have the same variance
    parametric = validate_double_parametric(min_n1, dpr1, min_n2, dpr2, alpha, log)
    
    dpr = dpr1 + dpr2
    test_group(dpr, alpha, parametric, log)
    test_groups(dpr1, dpr2, alpha, parametric, log)

def differences_across_classes(log, data, section, replicates, pooled_repeats=False, alpha=0.05):
    """
    Compare the ROIs in one class with those in the other classes.
    Example:
        From the ROIS in the isthmus, the ROIs classed as 'narrow-end' 
        can be compared agains the ROIs classed as 'wide_lumen', 'wide-end', and 'narrow-lumen'
    """
    l = f"{'=='*20}\nComparing ROIs across classes in {section} {'POOLING REPEATS' if pooled_repeats else ''}\n"
    print(l)
    print('H0: the samples were drawn from a population with the same distribution')
    log.write(l)
    log.write('H0: the samples/group of samples were drawn from a population with the same distribution\n')

    pairs_of_classes = [p for p in combinations([t[0] for t in classes2['ordering']], 2)]

    for roi_class_1, roi_class_2 in pairs_of_classes:
        rois_class1 = [roi for roi in classes2[section][roi_class_1]]
        rois_class2 = [roi for roi in classes2[section][roi_class_2]]
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
            
            process_repeat(data, rep, pooled_repeats, section, section, 
                replicates, replicates, rois_class1, rois_class2, log, alpha)

# comparisons for two cross sections
def differences_two_sections(log, data, sections, versions, replicates, paired_classes=True, pooled_repeats=False, alpha=0.05):
    """
    Compares the ROIs of a given class in the first section against 
    ROIs of non-matching classes in the second section.
    """
    mc = 'in MATCHING classes in' if paired_classes else 'in NON-matching classes in'
    pr = 'POOLING REPEATS' if pooled_repeats else 'PER REPEAT'
    l = f"{'=='*20}\nComparing ROIs {mc} {sections} {pr}\n"
    print(l)
    print('H0: both samples were drawn from a population with the same distribution')
    log.write(l)
    log.write('H0: both samples/group of samples were drawn from a population with the same distribution\n')
    
    class_reference = ['','']
    if type(sections) != list:
        sections = [sections, sections]
        replicates = [replicates, replicates]
    else:
        log.write(f'Renaming sections {sections}\n')
        class_reference[0] = sections[0][:-1] if sections[0][-1] == '1' else sections[0]
        class_reference[1] = sections[1][:-1] if sections[1][-1] == '2' else sections[1]
        log.write(f'New names {class_reference}\n')
    if paired_classes:
        pairs_of_classes = [(t[0], t[0]) for t in classes2['ordering']]
    else:
        pairs_of_classes = [p for p in combinations([t[0] for t in classes2['ordering']], 2)]

    for roi_class_1, roi_class_2 in pairs_of_classes:
        rois_class1 = [roi for roi in classes2[class_reference[0]][roi_class_1]]
        rois_class2 = [roi for roi in classes2[class_reference[1]][roi_class_2]]
        log.write(f"{'=='*10}\n")
        log.write(f'Rois in {sections[0]} {roi_class_1} {rois_class1}\n')
        log.write(f'Rois in {sections[1]} {roi_class_2} {rois_class2}\n')
        print(f'Rois in {sections[0]} {roi_class_1} {rois_class1}')
        print(f'Rois in {sections[1]} {roi_class_2} {rois_class2}')
        
        for rep in range(1 if pooled_repeats else min(replicates)):
            if not pooled_repeats:
                print(f'------------ repeat {rep+1} ------------')
                log.write(f'------------ repeat {rep+1} ------------\n')

            process_repeat(data, rep, pooled_repeats, sections[0], sections[1],
                replicates[0], replicates[1], rois_class1, rois_class2, log, alpha)

def header_single(log, section, version, path, iteration, rois, rep, pooled_repeats,
                number_concentration, size_distribution, bins, column):
    print('Comparing single section:')
    print(f'   section: {section}')
    print(f'   version: {version}')
    print(f' iteration: {iteration}')
    print(f'      rois: {rois}')
    print(f'replicates: {rep}')
    print(f'      path: {path}')
    print('------------------------------------')
    print('| pooled repeats?', pooled_repeats)
    print('| plot number concentration?', number_concentration)
    print('| plot size distribution?', size_distribution)
    print('| number of bins:', bins, 'number of columns:', columns)
    print('------------------------------------')
    log.write('Comparing single section:\n')
    log.write(f'   section: {section}\n')
    log.write(f'   version: {version}\n')
    log.write(f' iteration: {iteration}\n')
    log.write(f'      rois: {rois}\n')
    log.write(f'replicates: {rep}\n')
    log.write(f'      path: {path}\n')
    log.write('------------------------------------\n')
    log.write(f'| pooled repeats? {pooled_repeats}\n')
    log.write(f'| plot number concentration? {number_concentration}\n')
    log.write(f'| plot size distribution? {size_distribution}\n')
    log.write(f'| number of bins: {bins}, number of columns: {columns}\n')
    log.write('------------------------------------\n')

def header_double(log, section1, version1, path1, iter1, rois1, rep1, 
                section2, version2, path2, iter2, rois2, rep2, paired_classes, pooled_repeats, 
                number_concentration, size_distribution, bins, columns):
    print( 'Comparing->        [1]   |   [2]')
    print(f'   section:      {section1} | {section2}')
    print(f'   version:  {version1} | {version2}')
    print(f' iteration:        {int(iter1)/1e6:>.1f}M | {int(iter2)/1e6:>.1f}M')
    print(f'      rois:           {rois1} | {rois2}')
    print(f'replicates:            {rep1} | {rep2}')
    print(f'path [1]: {path1}')
    print(f'path [2]: {path2}')
    print('------------------------------------')
    print('| paired classes?', paired_classes)
    print('| pooled repeats?', pooled_repeats)
    print('| plot number concentration?', number_concentration)
    print('| plot size distribution?', size_distribution)
    print('| number of bins:', bins, 'number of columns:', columns)
    print('------------------------------------')
    log.write( 'Comparing->        [1]   |   [2]\n')
    log.write(f'   section:      {section1} | {section2}\n')
    log.write(f'   version:  {version1} | {version2}\n')
    log.write(f' iteration:        {int(iter1)/1e6:>.1f}M | {int(iter2)/1e6:>.1f}M\n')
    log.write(f'      rois:           {rois1} | {rois2}\n')
    log.write(f'replicates:            {rep1} | {rep2}\n')
    log.write(f'path [1]: {path1}\n')
    log.write(f'path [2]: {path2}\n')
    log.write('------------------------------------\n')
    log.write(f'| paired classes? {paired_classes}\n')
    log.write(f'| pooled repeats? {pooled_repeats}\n')
    log.write(f'| plot number concentration? {number_concentration}\n')
    log.write(f'| plot size distribution? {size_distribution}\n')
    log.write(f'| number of bins: {bins}, number of columns: {columns}\n')
    log.write('------------------------------------\n')


def compare_single(log, section, version, path, iteration, rois, rep, pooled_repeats,
                number_concentration, size_distribution, bins, column, across=False):
    """
    Compares the ROIs in each class wtih each other.
    Multiple groups are compared with a single test, if significant differences
    are found, pair-wise comparisons are performed
    """
    header_single(log, section, version, path, iteration, rois, rep, pooled_repeats,
                number_concentration, size_distribution, bins, column)

    source = f'./resources/analysis/output/stats_{version}_{iteration}_rois_v{rois}.pickle.bz2'
    print('Loading stats data from', source)
    with bz2.BZ2File(source, 'rb') as stats_source:
        s = pickle.load(stats_source)
        
    evs_per_replicate_df = s['evs_per_replicate_df']
    print(evs_per_replicate_df.head())

    #rep = 0
    d = evs_per_replicate_df.query('replicate==1')['radius_um']
    print(d[0])

    df = s['data_frame']
    
    c = df.query(f"section=='{section}' & replicate==1")['radius_um'].count()
    print(f'EVs in ROIS in replicate 1: {c}')
    log.write(f'EVs in ROIS in replicate 1: {c}\n')

    differences_within_classes(log, df, section, rep, pooled_repeats=pooled_repeats)
    
    if across:
        differences_across_classes(log, df, section, rep, pooled_repeats=pooled_repeats)

    # produce plots using the data in stats which comes with the dataframes
    # number concentration per ROI
    if number_concentration:
        plotting.number_concentration_per_roi(s, section, version, iteration, classes2[section], classes2['ordering'])

    if size_distribution:
        plotting.size_distribution_per_roi_np(s, section, version, iteration,
            classes2[section], rep, classes2['ordering'], columns=columns)

        from ev_model import plotting_mpl
        plotting_mpl.size_distribution_per_section(evs_per_replicate_df, section, version, iteration, rep)
        

def compare_double(log, section1, version1, path1, iter1, rois1, rep1, 
                section2, version2, path2, iter2, rois2, rep2, paired_classes, pooled_repeats, 
                number_concentration, size_distribution, bins, columns):
    
    header_double(log, section1, version1, path1, iter1, rois1, rep1, 
                section2, version2, path2, iter2, rois2, rep2, paired_classes, pooled_repeats, 
                number_concentration, size_distribution, bins, columns)
    stats = []
    data = None
    source1 = f'./resources/analysis/output/stats_{version1}_{iter1}_rois_v{rois1}.pickle.bz2'
    source2 = f'./resources/analysis/output/stats_{version2}_{iter2}_rois_v{rois2}.pickle.bz2'
    for source in [source1, source2]:
        with bz2.BZ2File(source, 'rb') as stats_source:
            s = pickle.load(stats_source)

            df = s['data_frame']
            
            if data is None:
                if section1 == section2:
                    print('adding 1 to the section name')
                    df.index = df.index.set_levels(df.index.levels[0].str.cat(['1']), level=0)
                data = df
            else:
                if section1 == section2:
                    print('adding 2 to the section name')
                    df.index = df.index.set_levels(df.index.levels[0].str.cat(['2']), level=0)
                data = data.append(s['data_frame'])
            stats.append(s)
    if section1 == section2:
        s1 = section1 + '1'
        s2 = section2 + '2'
    else:
        s1 = section1
        s2 = section2
    
    c1 = data.query(f"section=='{s1}' & replicate==1")['radius_um'].count()
    c2 = data.query(f"section=='{s2}' & replicate==1")['radius_um'].count()
    print(f'EVs in ROIS in replicate 1 from [1]: {c1}')
    print(f'EVs in ROIS in replicate 1 from [2]: {c2}')
    log.write(f'EVs in ROIS in replicate 1 from [1]: {c1}\n')
    log.write(f'EVs in ROIS in replicate 1 from [2]: {c2}\n')

    differences_two_sections(log, data, [s1, s2], [version1, version2], [rep1, rep2], 
            paired_classes=paired_classes, pooled_repeats=pooled_repeats)

    # produce plots using the data in stats which comes with the dataframes
    # number concentration per ROI
    if number_concentration:
        plotting.number_concentration_per_roi(stats[0], section1, version1, iter1, classes2[section1], classes2['ordering'])
        plotting.number_concentration_per_roi(stats[1], section2, version2, iter2, classes2[section2], classes2['ordering'])

    # size distribution per ROI
    if size_distribution:
        plotting.size_distribution_per_roi_np(stats[0], section1, version1, iter1,
            classes2[section1], rep1, classes2['ordering'], bins=bins, columns=columns)
        plotting.size_distribution_per_roi_np(stats[1], section2, version2, iter2,
            classes2[section2], rep2, classes2['ordering'], bins=bins, columns=columns)

if __name__ == '__main__':
    import argparse
    from datetime import datetime
    """
    comparisons possible:
    1 section - 1) within classes, 2) across classes
    2 sections - 1) within classes, 2) across classes
    """

    parser = argparse.ArgumentParser(add_help=True)
    s1 = parser.add_argument_group('section1', 'Settings for the first section to compare')
    s1.add_argument('section', choices=sections, help='The 1st section to analyse.')
    s1.add_argument('version', help="Input data version to read. Forms part of the file name to load: 'stats_SECTION_VERSION_ROIS.pickle.bz2 ")
    s1.add_argument('path', metavar='base_path', help='Directory where the target files are located (1st section to compare)')
    s1.add_argument('--rois', help='Version of the ROIs used for producing the stats')
    s1.add_argument('--iteration1', nargs='?', type=int, help='Iteration number to read - 1st section. Default 14.4k')
    s1.add_argument('--replicates1', nargs='?', type=int, help='The number of replicates to read - 1st section. Default 3')
    s1.add_argument('--pooled', action='store_true', help='Enables pooling the data from the repeats')
    s1.add_argument('--plotNC', action='store_true', help='Enables plotting the number concentration per ROI')
    s1.add_argument('--plotSD', action='store_true', help='Enables plotting the size distribution per ROI')
    s1.add_argument('--columns', nargs='?', type=int, help='Specify the number of columns of plots. Default: 3')
    s1.add_argument('--bins', nargs='?', type=int, help='Specify the number of bins for the histograms. Default: 15')
    
    s2 = parser.add_argument_group('section2', 
        """Settings for a second section to compare.
        Section2 is optional, however, it's required arguments are:
        --other_section, --other_roi_class, --other_path
        """)
    s2.add_argument('--section2', choices=sections, help='The 2nd section to analyse.')
    s2.add_argument('--version2', help="Input data version to read. Forms part of the file name to load: 'stats_SECTION_VERSION.pickle.bz2 ")
    s2.add_argument('--rois2', help='Version of the ROIs used for producing the stats (2nd section to compare)')
    s2.add_argument('--path2', help='Directory where the target files are located (2nd section to compare)') 
    s2.add_argument('--iteration2', nargs='?', type=int, help='Iteration number to read - 2nd section. Default 14.4k')
    s2.add_argument('--replicates2', nargs='?', type=int, help='The number of replicates to read - 2nd section. Default 3')
    s2.add_argument('--pairedClasses', action='store_true', help='If present, only ROIs from both sections in matching classes are compared')
    args = parser.parse_args()

    iter1 = args.iteration1 if args.iteration1 else 14400000
    i1 = f'{iter1/1e6:.1f}'
    r1 = args.rois
    rep1 = args.replicates1 if args.replicates1 else 3
    columns = args.columns if args.columns else 3
    bins = args.bins if args.bins else 15
    
    tstamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    p = '-pooled' if args.pooled else ''
    s, v = args.section, args.version

    if args.section2 or args.path2:
        iter2 = args.iteration2 if args.iteration2 else 14400000
        i2 = f'{iter2/1e6:.1f}'
        r2 = args.rois2
        rep2 = args.replicates2 if args.replicates2 else 3
        if not args.version2:
            raise RuntimeError('version2 not provided')
        
        pc = '-pairedClasses' if args.pairedClasses else '-nonPairedClasses'
        s2 = args.section2
        v2 = args.version2
        target_log = f"./resources/analysis/output/log_analysis_[{v}_{i1}M_{r1}-{v2}_{i2}M_{r2}]{p}{pc}_{tstamp}.txt"
        with open(target_log, 'w') as log:
            compare_double(log, args.section, args.version, args.path, iter1, r1, rep1,
                args.section2, args.version2, args.path2, iter2, r2, rep2, 
                args.pairedClasses, args.pooled, args.plotNC, args.plotSD, bins, columns)
    else:
        target_log = f"./resources/analysis/output/log_analysis_[{v}_{i1}M_{r1}]{p}_{tstamp}.txt"
        with open(target_log, 'w') as log:
            compare_single(log, args.section, args.version, args.path, iter1, r1, rep1, 
                    args.pooled, args.plotNC, args.plotSD, bins, columns)
    print('Log saved to:', target_log)