import microplates.plot as mp_plot
from IPython.display import display
import matplotlib.pyplot as plt

def get_samples_to_repeat(samples, gel,
                          repeat_if_status=['bad','?'], repeat_plates=None,
                          override_col = 'Override',
                          merge_on=['Plate','Well'],
                          require_cols=['Expt','Round','Sample']):

    bad_gel_samples = gel.loc[gel['Status'].isin(repeat_if_status),:]
    if override_col is not None:
        bad_gel_samples = bad_gel_samples[bad_gel_samples[override_col].str.strip() == '']
    if repeat_plates is not None:
        bad_gel_samples = bad_gel_samples.loc[bad_gel_samples['Plate'].isin(repeat_plates)]

    bad_samples = (pd.merge(left=bad_gel_samples, right=samples, how='left', on=merge_on)
                   .drop_duplicates(merge_on))
    bad_samples = bad_samples[~bad_samples[require_cols].isnull().any(axis='columns')]

    return bad_samples

def summarize_bad_samples(bad_samples):
    display(bad_samples['Status'].value_counts().to_frame())
    display(bad_samples['Phage Library'].value_counts().astype(int).to_frame())
    display(bad_samples['i5 bc group'].value_counts().to_frame())
    display(bad_samples['i7 bc group'].value_counts().to_frame())
    mp_plot.plate_heatmap(bad_samples['Well'].value_counts().rename('count').to_frame().reset_index().rename(columns={'index':'well'}),
                          parameter='count', annot=True)
    plt.figure()
    bad_samples['Plate'].plot.hist()

summarize_bad_samples(bad_samples)

from natsort import order_by_index, index_natsorted

def plate_repeat_samples(bad_samples,
                         libraries = ['Alpaca', 'Synthetic'],
                         start_wells = None,
                         separate_plates_per_library = True,
                         sort_by_columns = ['Phage Library','Expt','Round','Sample'], wells=96):

    repeat_samples = bad_samples #bad_samples[['Expt','Round','Sample','Phage Library','Status']]
    repeat_samples = repeat_samples.reindex(
        index=order_by_index(repeat_samples.index,
                             index_natsorted(zip(*(repeat_samples[x] for x in sort_by_columns)))))


    repeats_by_library = [repeat_samples[repeat_samples['Phage Library'] == x].copy() for x in libraries]
    if separate_plates_per_library:
        if start_wells is None: start_wells = ['A1' for r in libraries]
        for i,r in enumerate(repeats_by_library):
            r['Well'] = list(plates.utils.iterate_wells(len(r),    start = start_wells[i]))
    else:
        start_well = 'A1' if start_wells is None else start_wells
        for i,r in enumerate(repeats_by_library):
            r['Well'] = list(plates.utils.iterate_wells(len(r),    start = start_well))
            start_well = next(plates.utils.iterate_wells(1, start_well))
    return repeats_by_library

def print_repeats_table(repeat_alpaca, repeat_synthetic):
    with pd.option_context("display.max_rows", None, "display.max_columns", None):
        display(repeat_alpaca)
        display(repeat_synthetic)


def plot_repeat_plate(repeats, title):
    plt.figure(figsize=(12,8))
    mp_plot.plot_plate(repeats.set_index('Well'),
                           labels=['Expt','Round','Sample'], text_hue=['Expt','Round'], size=0)
    plt.suptitle(title);

def plot_cherrypick_bad_plates(bad_samples):
    for plate in sorted(bad_samples['Plate'].unique()):
        plt.figure()
        mp_plot.plot_cherrypick(bad_samples.loc[bad_samples['Plate'] == plate,'Well'])
        plt.suptitle("Plate {:.0f}".format(plate), y=1.05)
        plt.tight_layout()
