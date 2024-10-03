import microplates as plates
import microplates.plot as mp_plot
from IPython.display import display
import matplotlib.pyplot as plt
import pandas as pd
from natsort import order_by_index, index_natsorted

def add_analytical_gel_plate(date, plate, nrows=4, df=None, add_empty_row=True):
    rows = []
    for r in range(nrows):
        r1 = list(plates.iterwells(n=12, start=plates.row2letters(r*2)  +'1'))
        r2 = list(plates.iterwells(n=12, start=plates.row2letters(r*2+1)+'1'))
        
        rows.append({'Gel Date':date, 'Row':r+1, 'Lane':1, 'Plate':'', 'Notes':'Ladder'})
        i = 1
        for pair in zip(r1,r2):
             for well in pair:
                rows.append({'Gel Date':date, 'Row':r+1, 'Lane':i+1, 'Plate':plate, 'Well':well})
                i+=1
    if df is None:
        headers = "Gel Date	Row	Lane	Plate	Well	Status	Notes	Override	Override_Notes".split("\t")
        df = pd.DataFrame(columns=headers)
    df = pd.concat([df,pd.DataFrame(rows)])
    if add_empty_row:
        df = pd.concat([df, pd.DataFrame({'Gel Date':''},index=[0])])
    return df.fillna('')


def make_analytical_gel_template(plate_defs):
    """make a spreadsheet to track results of your analytical gel

    Parameters
    ----------
    plate_defs : list of dicts
        each dict represents one plate; each dict should have keys "date" and "plate"

    Example
    -------

    >>> df = make_analytical_gel_template([
    >>>     {'date':'2024-10-03', 'plate':1},
    >>>     {'date':'2024-10-03', 'plate':2},
    >>>     # etc.
    >>> ])
    >>> df.to_excel('./analytical_gel_template.xlsx', index=False)

    Returns
    -------
    pd.DataFrame
        spreadsheet template, with columns "Gel Date", "Row", "Lane", "Plate", "Well", "Status", "Notes", "Override", "Override_Notes".
    """
    df = None
    for plate in plate_defs:
        df = add_analytical_gel_plate(df=df, **plate)
    return df


def get_samples_to_repeat(samples, 
                          gel, 
                          repeat_if_status_is=['bad','?'], 
                          repeat_only_plates=None, 
                          override_col = 'Override',
                          merge_on=['Plate','Well'],
                          require_cols=['Expt','Round','Sample']):
    """figure out which samples need to have PCRs 1--2 repeated, based on results of analytical gel

    Parameters
    ----------
    samples : pd.DataFrame
        DataFrame of samples; should at least contain columns ['Plate','Well']
    gel : pd.DataFrame
        DataFrame containing results of analytical gel; should at least contain columns ['Plate','Well', 'Status, override_col]
    repeat_if_status_is : list, optional
        which values in the "Status" column should trigger samples to be repeated, by default ['bad','?']
    repeat_only_plates : list of int, optional
        if given, only repeat samples which appear in these library preparation plates, by default None
    override_col : str, optional
        if not None, any rows with a non-empty value in this column will NOT be repeated; by default 'Override'
    merge_on : list, optional
        which columns to merge on, by default ['Plate','Well']
    require_cols : list, optional
        drop any rows where any of these columns are NA, by default ['Expt','Round','Sample']

    Returns
    -------
    (pd.DataFrame, pd.DataFrame)
        tuple of two DataFrames: all_bad_samples, repeat_samples;
        - all_bad_samples : all samples that have a bad "Status" (i.e. Status is in `repeat_if_status_is`)
        - repeat_samples : all samples that should be repeated, i.e. they have a bad status, they are in `repeat_only_plates`, and they were not overridden
    """

    all_bad_gel_samples = gel.loc[gel['Status'].isin(repeat_if_status_is),:]
    all_bad_samples = (pd.merge(left=all_bad_gel_samples, right=samples, how='left', on=merge_on)
                   .drop_duplicates(merge_on))
    all_bad_samples = all_bad_samples[~all_bad_samples[require_cols].isnull().any(axis='columns')]
    
    repeat_samples = all_bad_samples
    if override_col is not None:
        repeat_samples = repeat_samples[repeat_samples[override_col].str.strip() == '']
    if repeat_only_plates is not None:
        repeat_samples = repeat_samples.loc[repeat_samples['Plate'].isin(repeat_only_plates)]
    
    return all_bad_samples, repeat_samples


def summarize_samples_to_repeat(repeat_samples):
    """summarize where the repeat samples came from

    Parameters
    ----------
    repeat_samples : pd.DataFrame
        result of calling `get_samples_to_repeat(...)[1]`
    """

    if len(repeat_samples) == 0:
        print("No samples to repeat!")
    
    display(repeat_samples['Status'].value_counts().to_frame())
    display(repeat_samples['Phage Library'].value_counts().to_frame())
    display(repeat_samples['Round'].value_counts().to_frame())
    display(repeat_samples['i5 bc group'].value_counts().to_frame())
    display(repeat_samples['i7 bc group'].value_counts().to_frame())
    plates.plot.plate_heatmap(repeat_samples['Well'].value_counts().rename('count').to_frame().reset_index().rename(columns={'index':'well'}),
                          parameter='count', annot=True)
    plt.xlabel('Column')
    plt.ylabel('Row')
    plt.title('# repeats per well position')
    plt.figure()
    repeat_samples['Plate'].plot.hist(title="# repeats per plate")


def plate_repeat_samples(repeat_samples, 
                         start_wells = None, 
                         separate_samples_by = 'Phage Library',
                         separate_into = 'plate',
                         sort_by_columns = ['Phage Library','Expt','Round','Sample'], 
                         wells=96):
    """arrange samples that are to be repeated on a plate or plates

    Use `separate_samples_by` and `separate_into` to arrange samples into separate plates or rows. 
    
    For example, in an experiment containing two values for `Phage Library`, `'Alpaca'` and `'Synthetic'`:

        >>> plate_repeat_samples(repeat_samples, separate_samples_by = 'Phage Library', separate_into = 'plates')
        # yields a list of DataFrames, one per plate, with one or more plates for 
        #'Alpaca' and one or more for 'Synthetic'.

        >>> plate_repeat_samples(repeat_samples, separate_samples_by = 'io', separate_into = 'plates')
        # yields a list of DataFrames, with input samples on a different plate 
        # from output samples.

        >>> plate_repeat_samples(repeat_samples, separate_samples_by = 'round', separate_into = 'rows')
        # yields a list of DataFrames, with samples from different rounds on different rows
    
    Parameters
    ----------
    repeat_samples : pd.DataFrame
        samples to repeat; should have columns [separate_samples_by] + sort_by_columns
    start_wells : list of str, optional
        for each plate of repeat samples, which well should they start on? if not given, each plate will start on well A1
    separate_samples_by : str, optional
        if given, separate samples which have different levels for this column, by default 'Phage Library'
    separate_into : str, optional
        how to separate levels of `separate_samples_by`: 
        - 'plate' (default) to place different levels on different plates; 
        - 'rows' to place each level on different rows
        - None to place samples on plate; do not separate by level of `separate_samples_by`. 
    sort_by_columns : list of str, optional
        sort the samples by these values, in order; by default ['Phage Library','Expt','Round','Sample']
    wells : int, optional
        number of wells in each plate, should be one onf 96, 384; by default 96

    Returns
    -------
    list of pd.DataFrame
        one DataFrame per plate, each comprised of rows of `bad_samples`; each DataFrame also has a column 'Well'.

    """

    repeat_samples = repeat_samples.reindex(
        index=order_by_index(repeat_samples.index, 
                             index_natsorted(zip(*(repeat_samples[x] for x in sort_by_columns)))))
    
    separate_repeats_by = separate_samples_by
    separate_repeats_by_values = repeat_samples[separate_repeats_by].unique()
    
    repeat_groups = [repeat_samples[repeat_samples[separate_repeats_by] == x].copy() for x in separate_repeats_by_values]
    
    if start_wells is None: start_wells = ['A1' for r in separate_repeats_by_values]
    if separate_into in ['plate', 'plates']:
        repeat_plates = []
        for i,r in enumerate(repeat_groups):
            total_wells = len(r)
            current_well = 0
            while (len(r) - current_well) > 0:
                last_well = min([current_well+wells, total_wells])
                n_wells = last_well-current_well
                
                r_plate = r.iloc[current_well:last_well].copy()
                r_plate['Well'] = list(plates.utils.iterate_wells(n_wells, start = start_wells[i]))
                repeat_plates.append(r_plate)
                
                current_well = current_well+wells
        return repeat_plates
    elif separate_into in ['row', 'rows']:
        start_well = start_wells[0]
        repeat_plates = []
        plate = []
        for i,r in enumerate(repeat_groups):
            _wells = list(plates.utils.iterate_wells(len(r),    start = start_well))
            r['Well'] = _wells
            last_well = plates.utils.cell2tuple(_wells[-1])
            plate.append(r)
            
            if last_well[0] >= plates.plate_layouts[wells][0]:
                repeat_plates.append(plate)
                plate = []
                start_well = 'A1'
            else:
                start_well = plates.utils.tuple2cell(last_well[0]+1, 0)
        if len(plate) > 0:
            repeat_plates.append(plate)
            
        return [pd.concat(plate) for plate in repeat_plates]
    else:
        repeat_plates = []
        r = repeat_samples
        total_wells = len(r)
        current_well = 0
        print(f"{total_wells} wells")
        while (len(r) - current_well) > 0:
            last_well = min([current_well+wells, total_wells])
            n_wells = last_well-current_well

            r_plate = r.iloc[current_well:last_well].copy()
            r_plate['Well'] = list(plates.utils.iterate_wells(n_wells))
            repeat_plates.append(r_plate)

            current_well = current_well+wells
            
        return repeat_plates

def print_repeats_table(*repeat_plates):
    with pd.option_context("display.max_rows", None, "display.max_columns", None):
        for plate in repeat_plates:
            display(plate)

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



def summarize_library(samples, 
                      normalization_volume = 5, 
                      normalization_mass = 12.5, 
                      reads_per_depth = 100_000, 
                      lane_cost = 2936, 
                      lane_reads = 2.8e9):
    """summarize the sequencing library

    Parameters
    ----------
    samples : pd.DataFrame
        samples included in the library
    normalization_volume : int, optional
        what volume, in microliters, will be drawn from each normalization well, by default 5 (assuming you add 10 microliters per well and use half)
    normalization_mass : float, optional
        what mass, in nanograms, will that normalization_volume contain, by default 12.5 (assuming each well recovers 20 nanograms and you use half)
    reads_per_depth : int, optional
        how many reads per sample, by default 100_000
    lane_cost : int, optional
        cost in dollars of one lane, by default $2936 (for NovaSeq X Plus 25B 2x150 flow cell at YCGA, as of 9/2024)
    lane_reads : float, optional
        how many reads per lane, by default 2.8e9 (for NovaSeq X Plus 25B 2x150 flow cell)

    Returns
    -------
    float, float
        number of reads, fraction of lane
    """
    mass = sum(samples['Depth'] * normalization_mass)
    volume = sum(samples.loc[samples['Depth'] == 1,'Depth'] * normalization_volume)
    depth = sum(samples['Depth'])
    reads = depth*reads_per_depth
    print(f"{len(samples)} total samples, with total depth {depth:g} = {reads:g} reads")
    print(f"{mass:>8g} ng total DNA")
    print(f"{volume:>8g} uL volume for 1-depth samples")
    for i,x in samples['Depth'].value_counts().iteritems():
        print(f" {i:>8g} depth : {x:>6g} samples")
    print(f"{reads / lane_reads:.3%} of a {lane_reads:g} lane; sequencing cost ~${lane_cost * reads / lane_reads:.2f}")
    return reads, (reads / lane_reads)