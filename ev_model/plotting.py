import holoviews as hv
import numpy as np
hv.extension('bokeh')
from holoviews import opts

def number_concentration_per_roi(d, section, version, iteration, rois_per_class, ordering):
    target = 'Total EVs'
    dt = []
    aggregateds = {}
    max_y = 0
    overlays = []
    for i in range(len(d['roi_areas'])):
        for j in range(len(d['total_evs_per_roi'][i])):
            # fetch the ROIs per class for this section
            for k,v in rois_per_class.items():
                if i+1 in v:
                        # class, roi, count
                    dt.append((k, i, d['total_evs_per_roi'][i][j]))
    dset = hv.Dataset(dt, ['Class', 'ROI', target], [target])
    for class_name, class_colour in ordering:
        aggregateds[class_name] = dset[class_name].aggregate(['ROI'], np.mean, np.std)
        if len(aggregateds[class_name][target]) > 0:
            m = max(aggregateds[class_name][target])
            if m > max_y:
                max_y = m
        overlays.append(hv.Bars(aggregateds[class_name]).opts(
            color=class_colour, fill_alpha=0.5, 
            title=class_name) * hv.ErrorBars(aggregateds[class_name]))
    
    # export this as an image
    bars = hv.Layout(overlays, 
            group=f'Number concentration per ROI, {section}', 
            label='').opts(opts.Bars(ylim=(0, max_y*1.1))).opts(toolbar=None)

    target_file = f"./resources/analysis/output/plot_{version}_{iteration}_numconc_per_roi.png"
    print('number_concentration_per_roi saved in:')
    print(target_file)
    hv.save(bars, target_file, fmt='png', size=100)

def size_distribution_per_roi_np(d, section, version, iteration, rois_per_class, replicates, 
        ordering, bins=10, columns=2):
    dt = []
    # create storage for the counts per roi and repeat
    for roi in range(len(d['roi_areas'])):
        l = []
        for rep in range(replicates):
            l.append(list())
        dt.append(l)
    
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
        fs = []    # storage for the frequencies per plot
        for rep in range(replicates):
            freqs, edges = np.histogram(dt[roi][rep], bins=bins, range=(mi,ma))
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
            if roi+1 in rois_per_class[cname]:
                col = color
                break
        h = hv.Histogram((edges, frequencies)).opts(color=col, fill_alpha=0.5, 
                title=f"ROI #{roi + 1} size dist.", xlabel='EV radius in um', 
                ylabel='Frequency', bgcolor="#E8DDCB", ylim=(0,maxy), padding=(0.05, 0))
        e = hv.ErrorBars((edges, frequencies, devs))
        plots.append((h * e))
    print('maxy:', maxy)
    target_file = f"./resources/analysis/output/plot_{version}_{iteration}_hist_sizedist_rois.png"
    print('Size distribution per ROI saved in:')
    print(target_file)
    hv.save(hv.Layout(plots).cols(columns).opts(toolbar=None), target_file, fmt='png', size=200)

    # exporting directly from bokeh works but depends on selenium and pillow
    # plus if using jupyter-lab: prior to launching jupyter-lab, execute: export OPENSSL_CONF=/etc/ssl/)
    # !conda install -c bokeh selenium -y
    # !conda install selenium pillow -y
    # !npm install -g phantomjs-prebuilt
    #
    #if produce_svg:
    #    render =  hv.render(hist_w_errors, backend='bokeh')
    #    render.output_backend = "svg"
    #    als.export_svgs(render, filename=f"./resources/analysis/output/{section_name}_{version}_roi_{i+1}.svg")

