import string
import random
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


def subplots(nrows=1, ncols=1, *, sharex=False, sharey=False, squeeze=True,
             width_ratios=None, height_ratios=None,
             subplot_kw=None, gridspec_kw=None, **fig_kw):
    from matplotlib.figure import Figure
    fig = Figure(**fig_kw)
    axs = fig.subplots(nrows=nrows, ncols=ncols, sharex=sharex, sharey=sharey,
                       squeeze=squeeze, subplot_kw=subplot_kw,
                       gridspec_kw=gridspec_kw, height_ratios=height_ratios,
                       width_ratios=width_ratios)
    return fig, axs

def display_all(*args):
    from IPython.display import display
    for a in args:
        display(a)

DEFAULT_REPLACEMENTS = {
    '∆{mexAB,oprM} ∆{mexCD,oprJ} ∆mexJKL ∆{mexHI,opmD} ∆opmH': '∆efflux1',
    '∆{mexAB,oprM} +nfxB ∆{mexCD,oprJ} ∆mexJKL ∆mexXY ∆opmH362 ∆{mexEF,oprN}': '∆efflux',
    'gacA::Tn(GmR) fliC::TcR ∆pilA': 'gacA::Tn',
    'retS::Tn(GmR) fliC::TcR ∆pilA': 'retS::Tn',
    '∆exsD attB::PexoT::ZTP_lacZ': 'ZTP riboswitch',
    'PAO1 PA0397': 'PA0397',
}

def shorten_descriptions(obs, replacements={}, new_column='desc_short', old_column='description'):
    from ..utils import replace_multiple
    replacements = {
        **DEFAULT_REPLACEMENTS,
        **replacements
    }
    obs[new_column] = replace_multiple(obs[old_column], replacements)
    return obs


def stars(pvalue):
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
    return "ns"

def trunc_ellipsis(max_len=50, ellipsis='...', split=None):
    if split is not None:
        def trunc(r):
            return split.join(x[:max_len]+ellipsis if len(x) > max_len else x for x in r.split(split))
    else:
        def trunc(x):
            return x[:max_len]+ellipsis if len(x) > max_len else x
    return trunc


def contrasting_color(hex_color):
    import colorsys
    from matplotlib import colors as mpl_colors
    hls = colorsys.rgb_to_hls(*mpl_colors.hex2color(hex_color))
    return '#000000' if hls[1] > 0.6 else '#ffffff'


def hash_to_color(h, **kwargs):
    others = {'others': '#aaaaaa', **kwargs}
    return others[h] if h in others else f'#{h[0:6]}'


def hash_to_mnemonicode(s, word_separator=' ', group_separator=',', format=True):
    import mnemonicode
    if s != 'others':
        if format:
            return mnemonicode.mnformat(bytes.fromhex(s), word_separator=word_separator, group_separator=group_separator)
        else:
            return list(mnemonicode.mnencode(bytes.fromhex(s)))
    else:
        return 'others'
hash_to_mn = hash_to_mnemonicode


def hash_to_mn_short(s, word_separator=' '):
    return hash_to_mnemonicode(s[:8], word_separator=word_separator)


def mn_to_hash(mn, word_separator=' ', group_separator=','):
    import mnemonicode
    return mnemonicode.mndecode([tuple(g.split(word_separator)) for g in mn.split(group_separator)]).hex()


def hash_to_proquint(s):
    import proquint
    return proquint.hex2quint_str(s[:8])


def hash_to_phonobyte(s):
    import phonobyte
    return phonobyte.encode(bytes.fromhex(s))


def pretty_hex(hex_str, sep=' ', n=2, length=6):
    from textwrap import wrap
    return sep.join(wrap(hex_str[0:length], n))


def pack_hex(hex_str):
    x = hex_str.replace(' ', '')
    if x[0] != '#':
        x = '#'+x
    return x


def contrasting_color(hex_color):
    from matplotlib.colors import rgb_to_hsv, hex2color
    hsl = rgb_to_hsv(hex2color(hex_color))
    return '#000000' if hsl[2] > 0.5 else '#ffffff'


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def png_to_html(png):
    return f"<img src='data:image/png;base64, {png}' />"

def repr_html_or_text(obj):
    _repr = ''
    # from IPython.display import HTML
    if hasattr(obj,'_repr_mime_'):
        for mime in ['text/html', 'text/markdown', 'image/png', 'text/plain']:
            _repr = obj._repr_mime_(mime)
            if _repr is not None:
                if mime == 'image/png':
                    return png_to_html(_repr)
                return _repr
    if hasattr(obj, '_repr_html_'):
        _repr = obj._repr_html_()
        if isinstance(_repr, tuple):
            _repr = _repr[0]
        if _repr is not None:
            return _repr
    if hasattr(obj, '_repr_markdown_'):
        _repr = obj._repr_markdown_()
        if isinstance(_repr, tuple):
            _repr = _repr[0]
        if _repr is not None:
            return _repr
    if hasattr(obj, '_repr_png_'):
        _png = obj._repr_png_()
        if _png is not None:
            return png_to_html(_png)
    if isinstance(obj, str):
        return obj
    try:
        if len(obj) > 0:
            if hasattr(obj[0], '_repr_html_') or hasattr(obj[0], '_repr_markdown_'):
                outs = [repr_html_or_text(o) for o in obj]
                return "\n".join(str(o) for o in outs if o is not None)
    except(TypeError):
        pass
    else:
        return repr(obj)


def page_break():
    from IPython.display import HTML, display
    display(HTML('<p style="page-break-after:always;"></p>'))


def pre(x):
    from IPython.display import HTML
    return HTML(f"<pre style='white-space:pre;'>{x}</pre>")


def display_table(rows):
    from IPython.display import HTML, display
    out = "<table>"
    for row in rows:
        out += "<tr>"
        for cell in row:
            out += "<td>"
            _repr = cell._repr_html_()
            if isinstance(_repr, tuple):
                _repr = _repr[0]
            out += _repr
            out += "</td>"
        out += "</tr>"
    return HTML(out)


def display_accordion(obj, title):
    from IPython.display import HTML, display
    if callable(obj):
        out = obj()
    else:
        out = obj

#     <div style='border: solid 1px silver;'><h2>{title}</h2>
#     {obj._repr_html_()}
#     </div>

#     display(HTML(f"""
#     <div class="lm-Widget p-Widget lm-Panel p-Panel p-Accordion jupyter-widgets widget-accordion widget-container">
#         <div class="lm-Widget p-Widget p-Collapse p-Accordion-child p-Collapse-open p-Accordion-child-active">
#             <div class="lm-Widget p-Widget p-Collapse-header">
#                 <span>{title}</span>
#             </div>
#             <div class="lm-Widget p-Widget lm-Panel p-Panel p-Collapse-contents">
#                  {obj._repr_html_()}
#             </div>
#         </div>
#     </div>
#     """))

    _id = id_generator()

    obj_repr = repr_html_or_text(obj)
    title_repr = repr_html_or_text(title)
    return (HTML(f"""
        <span class="trust-warning">Trust notebook to view interactive content.</span>
        <div class="tabs">
          <div class="tab">
            <input type="checkbox" id="{_id}1" name="{_id}">
            <label class="tab-label" for="{_id}1">{title_repr}</label>
            <div class="tab-content">
             {obj_repr}
            </div>
          </div>
        </div>
    """))


def show_all_accordions():
    from IPython.display import HTML, display
    display(HTML("""<style type='text/css'>
    input:checked ~ .tab-content {
        /* max-height: 100vh; */
        max-height: 100%;
        padding: 0.5em;
    }
    </style>"""))


def setup_accordion():
    from IPython.display import HTML, display
    display(HTML("""<style type='text/css'>
    .trust-warning {
      display: none;
    }
    .tabs>.tab>input[type=checkbox] {
      position: absolute;
      opacity: 0;
      z-index: -1;
    }
    /* Accordion styles */
    .tabs {
        border: solid 1px silver;
        overflow: hidden;
    }
    .tab {
      width: 100%;
      overflow: hidden;
    }
    .tab-label {
        display: flex;
        justify-content: space-between;
        padding: 0.5em;
        font-weight: bold;
        cursor: pointer;
    }
    /* Icon */
    .tab-label:hover {
      background: silver;
    }
    .tab-label::after {
      content: "\\25C0";
      width: 1em;
      height: 1em;
      text-align: center;
      transition: all .35s;
    }
    .tab-content {
        max-height: 0;
        padding: 0 0.5em;
        transition: all .35s;
        overflow-y: hidden;
    }

    // :checked
    input:checked + .tab-label {
        background: gray;
    }
    input:checked + .tab-label::after {
      transform: rotate(-90deg);
    }

    input:checked ~ .tab-content {
        /* max-height: 100vh; */
        max-height: 100%;
        padding: 0.5em;
    }
    </style>"""))


def repr_captured_output(co):
    from IPython.display import display, HTML
    return HTML(
                    (f"<pre style='font-size:smaller;'>{co.stderr}</pre>" if len(co.stderr) > 0 else '') +
                    (f"<pre>{co.stdout}</pre>" if len(co.stdout) > 0 else '') +
                    repr_html_or_text(co.outputs))


from contextlib import contextmanager
@contextmanager
def capture_accordion(title, **kwargs):
    from IPython.utils.capture import capture_output
    from IPython.display import display, HTML
    try:
        with capture_output() as co:
            yield co
    finally:
        display(
            display_accordion(
                repr_captured_output(co), 
                title, **kwargs)
        )

def rich_table(df,
               highlight_aa=['CDR', 'CDR1', 'CDR2',
                             'CDR3', 'AA', 'AA_aligned'],
               highlight_na=['NA', 'sequence', 'aligned', 'NA_aligned'],
               monospace=[], color=['color'], auto_color=True, optimize=True,
               caption=None, format=None, **kwargs):
    from pandas.io.formats.style import Styler

    df = df.copy()
    if auto_color and len(color) > 0:
        color_col = ['CDR3ID', 'aaSVID', 'ASVID']

        if df.index.name in color_col:
            df['color'] = df.index.str.slice(0, 6)
        else:
            for col in color_col:
                if col in df.index.names:
                    df['color'] = df.index.get_level_values(
                        col).str.slice(0, 6)
                    break
                elif col in df.columns:
                    df['color'] = df[col].str.slice(0, 6)
                    break

    columns = set(df.columns)

    highlight_aa = set(highlight_aa) - set(highlight_na)
    aa_highlight_columns = columns.intersection(highlight_aa)
    if len(aa_highlight_columns) > 0:
        from .syntax import aa_highlighter
        for col in aa_highlight_columns:
            df[col] = df[col].apply(lambda x: aa_highlighter.highlight(x))

    na_highlight_columns = columns.intersection(highlight_na)
    if len(na_highlight_columns) > 0:
        from .syntax import na_highlighter
        for col in na_highlight_columns:
            df[col] = df[col].apply(lambda x: na_highlighter.highlight(x))

    color_columns = columns.intersection(color)
    monospace_columns = columns.intersection(
        set(monospace) | aa_highlight_columns | na_highlight_columns)

    if optimize:
        # uuid_len=0 will prevent the table from getting a UUID and make styles
        # conflict across multiple tables
        dfs = Styler(df, cell_ids=False)  # uuid_len=0, cell_ids=False)
    else:
        dfs = df.style
    
    if format is not None:
        dfs = dfs.format(format)

    table_styles = []
    column_index = dict(zip(list(df.columns), range(len(df.columns))))
    if len(color_columns) > 0:
        # dfs = dfs.applymap(lambda c: f"background-color: #{c}; color: {contrasting_color('#'+c)}", subset=color_columns)
        # nb: applymap deprecated in pandas ~2.1, replaced by map(axis=None)
        dfs = dfs.applymap(lambda c: f"background-color: #{c}; color:#ffffff",
                           # subset must be a valid input to DataFrame.loc if more than one column is given,
                           # hence this is basically dfs.loc[:,color_columns]
                           # can't index using a set, so cast to list
                           subset=(slice(None), list(color_columns)))
    if len(monospace_columns) > 0:
        for col in monospace_columns:
            table_styles.append(
                {'selector': f'td.col{column_index[col]}, td.col{column_index[col]} pre', 'props': "font-family:monospace;white-space:pre;"})

        # dfs = dfs.applymap(lambda c: f"font-family:monospace;white-space:pre;",
        #                    subset=(slice(None), monospace_columns))
    dfs = dfs.set_table_styles(table_styles)

    if caption is not None:
        dfs = dfs.set_caption(caption)
        
    return dfs


def tag(id, space='cdr3', library='??', mn=True, short=False, html=True):
    mns = hash_to_mn(id, format=False) #list(mnemonicode.mnencode(bytes.fromhex('7fd62864f5fa6e2cc2076c4c046eb109')))

    bg = hash_to_color(id)
    fg = contrasting_color(bg)
    if library is None:
        library = '??'
    if short:
        out = (
            "<span>"
            f"<small>{space.upper()}</small>&nbsp;"
            f"<span style='background-color:{bg};color:{fg};padding:0.3em;'>{pretty_hex(id)}</span>&nbsp;"
            f"\"<b>{' '.join(mns[0])}</b>\""
            "</span>"
        )
    else:
        out = (
            "<div style='margin:0.5em 0;border-bottom:solid 1px silver;'>"
            f"<small>{space.upper()}</small>&nbsp;"
            f"<span style='background-color:{bg};color:{fg};padding:0.3em;'>{pretty_hex(id)}</span>&nbsp;"
            f"<small style='font-size:8pt;'><code>{id}</code>&nbsp;({library})</small><br />"
            f"<b>{'_'.join(mns[0])}</b> {' '.join('_'.join(m) for m in mns[1:])}"
            "</div>"
        )
    if html:
        from IPython.display import display, HTML
        return HTML(out)
    else:
        return out



def cdr3_tag(CDR3ID, library=None):
    return tag(id=CDR3ID, space='cdr3', library=library)


def cleanup_patchwork():
    import matplotlib.pyplot as plt
    import patchworklib as pw
    plt.close(pw.Brick._figure)
    pw.Brick._figure = plt.figure(figsize=(1, 1))


def pprint_dict(d, number_format='{:.2g}'):
    def fmt(k, v):
        if isinstance(v, int):
            v = str(v)
        elif isinstance(v, float):
            v = number_format.format(v)
        elif isinstance(v, str):
            v = f"'{v}'"
        return f"'{k}': {v}"
    return '{' + ', '.join(fmt(k, v) for k, v in d.items()) + '}'


def histplot_sample(X, n=1000, **kwargs):
    import seaborn as sns
    import numpy as np
    return sns.histplot(data=np.random.choice(X, n), **kwargs)


def plot_grid(rows, cols, f, subplot_width=4, figsize=None, **kwargs):
    import matplotlib.pyplot as plt
    if not isinstance(subplot_width, tuple):
        subplot_width = (subplot_width, subplot_width)
    if figsize is None:
        figsize = (subplot_width[0]*len(cols), subplot_width[1]*len(rows))
    fig, axs = plt.subplots(ncols=len(cols), nrows=len(rows),
                            figsize=figsize,
                            layout='tight', **kwargs)
    for j, col in enumerate(cols):
        for i, row in enumerate(rows):
            f(row, col, ax=axs[i, j])


def css_bar(start: float, end: float, color: str) -> str:
    """
    Generate CSS code to draw a bar from start to end in a table cell.
    Uses linear-gradient.

    from https://github.com/pandas-dev/pandas/blob/main/pandas/io/formats/style.py#L3786

    Parameters
    ----------
    start : float
        Relative positional start of bar coloring in [0,1]
    end : float
        Relative positional end of the bar coloring in [0,1]
    color : str
        CSS valid color to apply.
    Returns
    -------
    str : The CSS applicable to the cell.
    Notes
    -----
    Uses ``base_css`` from outer scope.
    """
    cell_css = ''
    if end > start:
        cell_css += "background: linear-gradient(90deg,"
        if start > 0:
            cell_css += f" transparent {start*100:.1f}%, {color} {start*100:.1f}%,"
        cell_css += f" {color} {end*100:.1f}%, transparent {end*100:.1f}%)"
    return cell_css


def rug_bar(x, vmin, vmax, color='#d65f5f'):
    """generate CSS for a bar starting from zero, extending to x, where vmin < x < vmax

    Parameters
    ----------
    x : float
        position of bar, in [vmin, vmax]
    vmax : float
        max value
    vmin : float
        min value
    color : str, optional
        what color to make the bar; by default '#d65f5f'

    Returns
    -------
    str
        CSS string describing the bar
    """
    rng = vmax - vmin
    zero_point = abs(vmin) / rng
    if x > 0:
        start = zero_point
        end = (x-vmin)/rng
    elif x < 0:
        end = zero_point
        start = (x-vmin)/rng
    else:
        start = end = zero_point
    return css_bar(start, end, color)


def pickle_figure(fig, buf=None):
    """saves a matplotlib figure with pickle to the provided file handle or a BytesIO

    Parameters
    ----------
    fig : plt.Figure
    buf : file handle, optional
        if None given, `io.BytesIO` is created

    Returns
    -------
    file
        buffer to which the file was written
    """
    import pickle
    from io import BytesIO
    if buf is None:
        buf = BytesIO()
    pickle.dump(fig, buf)
    return buf


def figure_to_buffer(fig, buf=None, close=True):
    """writes a matplotlib figure to a file handle (or creates a BytesIO buffer), then closes figure and returns the file

    Parameters
    ----------
    fig : _type_
        _description_
    buf : file, optional
        if None given, `io.BytesIO` is created

    Returns
    -------
    file
    """
    import matplotlib.pyplot as plt

    buf = pickle_figure(fig, buf=buf)
    if close:
        plt.close(fig)
    else:
        fig.clf()
    return buf


def buffer_to_figure(buf):
    """loads a pickled matplotlib figure
    """

    import pickle
    buf.seek(0)
    fig = pickle.load(buf)
    return fig


def figure_to_sparkline(fig, axis_off=True, subplots_adjust=True, img_extra='', close=True):
    import matplotlib.pyplot as plt
    import base64
    from io import BytesIO

    if axis_off:
        for ax in fig.axes:
            ax.axis('off')

    # squeeze axis to the edges of the figure
    if subplots_adjust:
        fig.subplots_adjust(left=0)
        fig.subplots_adjust(right=0.99)
        fig.subplots_adjust(bottom=0.1)
        fig.subplots_adjust(top=0.9)

    # save the figure to html
    bio = BytesIO()
    fig.savefig(bio)
    html = f"""<img src="data:image/png;base64,{base64.b64encode(bio.getvalue()).decode('utf-8')}" {img_extra}/>"""
    if close:
        plt.close(fig)
    else:
        fig.clf()
    return html


def sparkline(data, plot, fig=None,
              figsize=(4, 0.25), sparkline_kwargs={}, **kwargs):
    """Create a single HTML image tag containing a base64 encoded
    sparkline style plot

    Parameters
    ----------
    data : array-like (list, 1d Numpy array, Pandas Series) sequence of
        data to plot
    figsize : tuple of float, length and height of sparkline plot.  Attribute
        of matplotlib.pyplot.plot.
    """
    import matplotlib.pyplot as plt

    # fig = plt.figure(figsize=figsize)  # set figure size to be small
    # ax = fig.add_subplot(111)
    if fig is not None:
        close=False
        # ax = fig.add_subplot(1,1,1)
        ax = fig.subplots()
        # ax = fig.axes[0]
    else:
        fig, ax = plt.subplots(figsize=figsize)
        close=True

    plot(data, ax=ax, fig=fig, **kwargs)

    return figure_to_sparkline(fig, close=close, **sparkline_kwargs)


def plot_mini_histogram(
    values, ax=None, fig=None, vmin=None, vmax=None, 
    marks=None, mark_kws={}, formatter=None, figsize=(2, 0.5), **kwargs):
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    import pandas as pd

    if ax is None:
        # fig = Figure(figsize=figsize)
        # ax = fig.subplots()
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()


    if marks is not None:
        for mark in marks:
            ax.axvline(mark, **mark_kws)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # ticks = [round(values.min(), 1), 0, round(values.max(), 1)]
    if vmax is None:
        vmax = round(values.max(), 1)
    if vmin is None:
        vmin = round(values.min(), 1)
    if vmax > 0 and vmin < 0:
        ticks = [vmin, 0, vmax]
        labels = [vmin, '', vmax]
    else:
        ticks = [vmin, vmax]
        labels = ticks

    if pd.isna(ticks[0]) and pd.isna(ticks[-1]):
        pass
    else:
        ax.hist(values, **kwargs)

    if formatter is not None:
        labels = [formatter.format(x) if x != '' else x for x in labels]

    ax.set_xticks(ticks, labels=labels)
    ax.set_yticks([])
    ax.tick_params(axis='x', direction='in', labelbottom=True, pad=-16)
    ax.margins(x=0.1)
    return fig


def mini_histogram(values, vmin=None, vmax=None, **kwargs):
    fig = plot_mini_histogram(values, vmin=vmin, vmax=vmax, **kwargs)
    return figure_to_sparkline(fig, axis_off=False, subplots_adjust=True)

# HTML(mini_histogram(df_enr.query(f"CDR3ID == '{feature}'")['log_enrichment']))
# mini_histogram(df_enr.query(f"CDR3ID == '{feature}'")['log_enrichment'])


def plot_abundance_sparkline(
        data, ax, fig,
        point=True, point_color='red', point_marker='.',
        point_fill='red', point_size=6, point_alpha=1.0,
        points=True, fill=True, fill_color='blue', fill_alpha=0.1, **kwargs):
    """
    point : bool, show point marker on last point on right
    point_location : not implemented, always plots rightmost
    point_color : str, matplotlib color code for point, default 'red'
    point_marker : str, matplotlib marker code for point
    point_fill : str, matplotlib marker fill color for point, default 'red'
    point_size : int, matplotlib markersize, default 6
    point_alpha : float, matplotlib alpha transparency for point
    fill : bool, show fill below line
    fill_color : str, matplotlib color code for fill
    fill_alpha : float, matplotlib alpha transparency for fill
    **kwargs : keyword arguments passed to matplotlib.pyplot.plot
    """
    import matplotlib.pyplot as plt

    data = list(data)

    plot_len = len(data)
    plot_min = min(data)
    point_x = plot_len - 1

    ax.plot(data, **kwargs)

    # fill between the axes
    ax.fill_between(range(plot_len), data, plot_len*[plot_min],
                     color=fill_color, alpha=fill_alpha)

    # plot the right-most point red, probably on makes sense in timeseries
    ax.plot(point_x, data[point_x], color=point_fill,
             marker=point_marker, markeredgecolor=point_color,
             markersize=point_size,
             alpha=point_alpha, clip_on=False)


def abundance_sparkline(data, **kwargs):
    return sparkline(data, plot=plot_abundance_sparkline, **kwargs)


def extract_encoded_data(chart, d=None, overwrite={}, extra_fields = [], view_name_col=None, verbose=False, datasets={}, data=None):
    """export data that appears in an altair chart as a Pandas DataFrame

    This is intended as an easy way to generate "source data" from an existing chart. Each data point will be represented by one row in the DataFrame. 
    Each field that is encoded to a channel in the chart (e.g. X, Y, color/fill, etc.) will be represented by a column in the DataFrame. If a "title" 
    is given for that encoding, that title will be the name of the column, otherwise it will be the original name in the source data.
    
    View composition operators 'facet', 'hconcat', and 'vconcat' are supported but 'repeat' is not.
    Transformations are not supported, nor are dynamic params (or selections). You can use `overwrite` to replace these computed fields with other columns 
    that are already in the dataset. For complex charts with multiple 'datasets' or Vega-Lite 'data' arguments, this function will do its best to find the 
    data source corresponding to a particular view, but it may not be perfect.

    If multiple views are composed (e.g. with hconcat or vconcat), one DataFrame will be generated for each view and the DataFrames will be vertically 
    concatenated. The columns in the final DataFrame will thus be the union of all columns in all views. You can use the `view_name_col` to add an 
    additional column listing the view title (e.g. set using `.properties(title='TITLE')`) or the auto-generated view 'name' (e.g. 'view_NN') from which each data
    point came. This is useful for charts with multiple subplots. You can could then separate out these DataFrames and drop extra columns like this:

        >>> df = extract_encoded_data(ch, view_name_col='subplot')
        >>> for title, dff in df.groupby('subplot'):
        ...     print(title)
        ...     display(dff.dropna(axis='columns', how='all'))

    To use this function for a Vega-lite plot within a dashboard visualization from nbseq.viz.dash, click the three dots and select "View Source," then copy
    the JSON source code into this snippet:

        >>> chart = alt.Chart.from_json('''JSON_SOURCE_HERE''')
        >>> extract_encoded_data(chart, view_name_col='subplot')

    Parameters
    ----------
    chart : alt.Chart
        chart produced by alt.Chart or loaded by alt.from_json, alt.from_dict, etc.
    d : dict, optional
        (private, used for recursion) the Vega-Lite chart specification as dict, or a subset thereof; if omitted, will be extracted from `chart`
    overwrite : dict, optional
        if given keys should be names of encoding fields that you would like to replace, and values
        should other columns in the data. For example, if color was mapped to an interactive parameter, 
        you can replace it with a column value; by default {}
    extra_fields : list, optional
        additional columns from the data to add to each view by default []
    view_name_col : str, optional
        if given, add an additional column with this name to the output; this column will contain the title or name of the view from which each data point came from

    Returns
    -------
    pd.DataFrame
        data
    """
    import pandas as pd
    import altair as alt
    import warnings

    def vprint(*args, **kwargs):
        if verbose: print(*args, **kwargs)
    
    
    if d is None:
        d = chart.to_dict()
    
    vprint(f"step: {str(d)[:80]}...")
    vprint(f"- data = {str(data)[:80]}")
    
    if 'datasets' in d:
        datasets = {
            **datasets, 
            **{name: pd.DataFrame(data) for name, data in d['datasets'].items()}
        }
        
    if 'data' in d:
        vprint(f"- d['data'] = {d['data']}")
        if 'name' in d['data']:
            if d['data']['name'] in datasets:
                data = datasets[d['data']['name']]
            else:
                raise ValueError(f"Cannot locate named dataset '{d['data']['name']}'")
        elif 'values' in d['data']:
            data = pd.DataFrame(d['data']['values'])
        else:
            raise NotImplementedError(f"Cannot process vega-lite data directive {d['data']}")
    if data is None:
        data = chart.data

    
    dfs = []
    for operator in ['vconcat', 'hconcat']:
        if operator in d:
            for c in d[operator]:
                dfs.append(extract_encoded_data(chart, c, overwrite=overwrite, extra_fields=extra_fields, view_name_col=view_name_col, verbose=verbose, datasets=datasets, data=data))
    if len(dfs) > 0:
        return pd.concat(dfs)
    
    fields = []
    
    def parse_encoding(o):
        channels = list(o.keys())

        # move these channels to the front of the line if present, that way those columns appear first in the data
        for channel in 'row', 'column', 'facet':
            if channel in channels: 
                channels.insert(0, channels.pop(channels.index(channel)))
                
        for channel in channels:
            encoding = o[channel]
            field = None
            title = None
            if 'field' in encoding:
                field = encoding['field']
            if 'title' in encoding and encoding['title'] is not None:
                title = encoding['title']
            else: title = field
                
            if field is not None:
                fields.append((field, title))
                
    if 'facet' in d:
        parse_encoding(d['facet'])
    if 'encoding' in d:
        parse_encoding(d['encoding'])
    if 'spec' in d:
        if 'encoding' in d['spec']:
            parse_encoding(d['spec']['encoding'])
        if 'layer' in d['spec']:
            for layer in d['spec']['layer']:
                if 'encoding' in layer:
                    parse_encoding(layer['encoding'])
    if 'title' in d:
        name = d['title']
    elif 'name' in d:
        name = d['name']
    else:
        name = None
        
    # if field name is in the "overwrite" dict, replace
    fields = [((overwrite[f], title) if f in overwrite else (f,title)) for f, title in fields]


    if data is alt.Undefined:
        if len(fields) > 0:
            raise ValueError(f"Could not locate data for view where fields {fields} are mapped; view (or parent views) did not have embedded 'data' statement and chart.data is Undefined")
        else:
            return pd.DataFrame()
    
    # add extra fields
    fields = fields + [(ef if isinstance(ef, tuple) else (ef, ef)) for ef in extra_fields]

    # deduplicate fields, preserving order
    fields = list(dict.fromkeys(fields))
    
    # filter fields to only those in the chart
    columns = []
    titles = []
    for (field, title) in fields:
        if field not in data:
            warnings.warn(f"Column '{field}' not found in chart data, removing")
        else:
            columns.append(field)
            titles.append(title)

    df = pd.DataFrame(data[columns])
    df.columns = titles
    vprint(f"- columns {columns}")
    vprint(f"- renamed {titles}")
    if (view_name_col is not None) and (name is not None):
        df.insert(0, view_name_col, name)
    vprint(f"return {str(df.columns)[:80]}...")
    return df