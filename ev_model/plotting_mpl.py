import holoviews as hlv
import numpy as np
# from holoviews import opts

hlv.extension('matplotlib', 'bokeh')


def size_distribution_per_section(data,
                                  section,
                                  version,
                                  iteration,
                                  replicates,
                                  bins=10):
    hlv.renderer('matplotlib')
    plots = []
    max_y, min_y = 0, None
    fs = []
    for rep in range(replicates):
        d = data.query(f'replicate=={rep}')[['radius_um']]
        freqs, edges = np.histogram(d, bins=bins)
        fs.append(freqs)
    # compute error bars per bin
    # compute mean, sd per bin
    print('edges:', edges)
    values = np.zeros([bins, replicates])
    for bi in range(bins):
        for rep in range(replicates):
            v = fs[rep][bi]
            if min_y is None or v < min_y:
                min_y = v
            max_y = v if v > max_y else max_y
            values[bi][rep] = fs[rep][bi]
    frequencies = values.mean(axis=1)
    devs = values.std(axis=1)
    # create the histograms
    h = hlv.Histogram((edges, frequencies)).opts(  # fill_alpha=0.5,
        title=f"{section} size dist.",
        xlabel='EV radius in um',
        ylabel='Frequency',
        bgcolor="#E8DDCB",
        padding=(0.02, 0),
        ylim=(0.9 * min_y, 1.1 * max_y),
        logy=True)  # , yticks=[1000,1500,3000,6000,90000,15000,18000])
    e = hlv.ErrorBars((edges, frequencies, devs))
    plots.append((h * e))

    print('max_y:', max_y, 'min_y:', min_y)
    target_file = f"./resources/analysis/output/plot_{version}_{iteration}_hist_sizedist_sections.png"
    print('Size distribution per section saved in:')
    print(target_file)
    hlv.save(hlv.Layout(plots), target_file, fmt='png', size=200)
