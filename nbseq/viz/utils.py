import string
import random
import warnings
warnings.filterwarnings("ignore",category=FutureWarning)


def trunc_ellipsis(max_len=50, ellipsis='...', split=None):
    if split is not None:
        def trunc(r):
            return split.join(x[:max_len]+ellipsis if len(x) > max_len else x for x in r.split(split))
    else:
        def trunc(x):
            return x[:max_len]+ellipsis if len(x) > max_len else x
    return trunc

def contrasting_color(hex_color):
    hls = colorsys.rgb_to_hls(*mpl_colors.hex2color(hex_color))
    return '#000000' if hls[1] > 0.6 else '#ffffff'

def hash_to_color(h, **kwargs):
    others = { 'others':'#aaaaaa', **kwargs }
    return others[h] if h in others else f'#{h[0:6]}'

def pretty_hex(hex_str, sep=' ',n=2, length=6):
    from textwrap import wrap
    return sep.join(wrap(hex_str[0:length], n))

def pack_hex(hex_str):
    x = hex_str.replace(' ','')
    if x[0] != '#':
        x = '#'+x
    return x


def contrasting_color(hex_color):
    from matplotlib.colors import rgb_to_hsv, hex2color
    hsl = rgb_to_hsv(hex2color(hex_color))
    return '#000000' if hsl[2] > 0.5 else '#ffffff'

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def repr_html_or_text(obj):
    _repr = ''
    if hasattr(obj,'_repr_html_'):
        _repr = obj._repr_html_()
        if isinstance(_repr, tuple):
            _repr = _repr[0]
    elif hasattr(obj,'_repr_markdown_'):
        _repr = obj._repr_markdown_()
        if isinstance(_repr, tuple):
            _repr = _repr[0]
    elif isinstance(obj, str):
        _repr = obj
    else:
        _repr = repr(obj)
    return _repr


def page_break():
    from IPython.display import HTML, display
    display(HTML('<p style="page-break-after:always;"></p>'))

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
        overflow: hidden;
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



def rich_table(df, highlight_aa=['CDR','CDR1','CDR2','CDR3','AA'], monospace=[], color=['color'], auto_color=True, optimize=True, **kwargs):
    from pandas.io.formats.style import Styler

    df = df.copy()
    if auto_color and len(color) > 0:
        color_col = ['CDR3ID','aaSVID','ASVID']

        if df.index.name in color_col:
            df['color'] = df.index.str.slice(0,6)
        else:
            for col in color_col:
                if col in df.columns:
                    df['color'] = df[col].str.slice(0,6)
                    break

    columns = set(df.columns)

    highlight_columns = columns.intersection(highlight_aa)
    if len(highlight_columns) > 0:
        from .syntax import aa_highlighter
        for col in highlight_columns:
            df[col] = df[col].apply(lambda x: aa_highlighter.highlight(x))

    color_columns = columns.intersection(color)
    monospace_columns = columns.intersection(set(monospace) | highlight_columns)

    if optimize:
        # uuid_len=0 will prevent the table from getting a UUID and make styles
        # conflict across multiple tables
        dfs = Styler(df, cell_ids=False) #uuid_len=0, cell_ids=False)
    else:
        dfs = df.style

    table_styles = []
    column_index = dict(zip(list(df.columns), range(len(df.columns))))
    if len(color_columns) > 0:
        # dfs = dfs.applymap(lambda c: f"background-color: #{c}; color: {contrasting_color('#'+c)}", subset=color_columns)
        dfs = dfs.applymap(lambda c: f"background-color: #{c}; color:#ffffff",
                           # subset must be a valid input to DataFrame.loc if more than one column is given,
                           # hence this is basically dfs.loc[:,color_columns]
                           subset=(slice(None),color_columns))
    if len(monospace_columns) > 0:
        for col in monospace_columns:
            table_styles.append({'selector': f'td.col{column_index[col]}, td.col{column_index[col]} pre', 'props':"font-family:monospace;white-space:pre;" })

        # dfs = dfs.applymap(lambda c: f"font-family:monospace;white-space:pre;",
        #                    subset=(slice(None), monospace_columns))
    dfs = dfs.set_table_styles(table_styles)
    return dfs


def tag(id, space='cdr3', library='??'):
    from IPython.display import display, HTML
    bg = hash_to_color(id)
    fg = contrasting_color(bg)
    if library is None:
        library = '??'
    display(HTML(f"<small>{space.upper()}</small>&nbsp;<span style='background-color:{bg};color:{fg};padding:0.3em;'>{pretty_hex(id)}</span>&nbsp;<small style='font-size:8pt;'><code>{id}</code>&nbsp;({library})</small>"))

def cdr3_tag(CDR3ID, library=None):
    return tag(id=CDR3ID, space='cdr3',library=library)

def cleanup_patchwork():
    import patchworklib as pw
    plt.close(pw.Brick._figure)
    pw.Brick._figure = plt.figure(figsize=(1,1))


def histplot_sample(X, n=1000, **kwargs):
    import seaborn as sns
    import numpy as np
    return sns.histplot(data=np.random.choice(X,n), **kwargs)

def plot_grid(rows, cols, f, subplot_width=4, figsize=None, **kwargs):
    if not isinstance(subplot_width, tuple):
        subplot_width = (subplot_width, subplot_width)
    if figsize is None:
        figsize=(subplot_width[0]*len(cols), subplot_width[1]*len(rows))
    fig, axs = plt.subplots(ncols=len(cols), nrows=len(rows),
                        figsize=figsize,
                        layout='tight', **kwargs)
    for j, col in enumerate(cols):
        for i, row in enumerate(rows):
            f(row, col, ax=axs[i,j])
