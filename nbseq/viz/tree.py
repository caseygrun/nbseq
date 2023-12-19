from .asv import hash_to_color, pretty_hex, pack_hex
from .utils import contrasting_color


def get_depths(tree, unit_branch_lengths=False):  # noqa: D402
    """Create a mapping of tree clades to depths (by branch length).
    :Parameters:
        unit_branch_lengths : bool
            If True, count only the number of branches (levels in the tree).
            By default the distance is the cumulative branch length leading
            to the clade.
    :returns: dict of {clade: depth}, where keys are all of the Clade
        instances in the tree, and values are the distance from the root to
        each clade (including terminals).
    """  # noqa: D402
    if unit_branch_lengths:
        depth_of = lambda c: 1  # noqa: E731
    else:
        depth_of = lambda c: c.length or 0  # noqa: E731
    depths = {}

    def update_depths(node, curr_depth):
        depths[node] = curr_depth
        for child in node.children:
            new_depth = curr_depth + depth_of(child)
            update_depths(child, new_depth)

    update_depths(tree, tree.length or 0)
    return depths

def draw(
    tree,
    label_func=str,
    do_show=True,
    show_confidence=True,
    # For power users
    axes=None,
    branch_labels=None,
    label_colors=None,
    *args,
    **kwargs,
):
    """Plot the given tree using matplotlib (or pylab).
    The graphic is a rooted tree, drawn with roughly the same algorithm as
    draw_ascii.
    Additional keyword arguments passed into this function are used as pyplot
    options. The input format should be in the form of:
    pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict), or
    pyplot_option_name=(dict).
    Example using the pyplot options 'axhspan' and 'axvline'::
        from Bio import Phylo, AlignIO
        from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
        constructor = DistanceTreeConstructor()
        aln = AlignIO.read(open('TreeConstruction/msa.phy'), 'phylip')
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        tree = constructor.upgma(dm)
        Phylo.draw(tree, axhspan=((0.25, 7.75), {'facecolor':'0.5'}),
        ... axvline={'x':0, 'ymin':0, 'ymax':1})
    Visual aspects of the plot can also be modified using pyplot's own functions
    and objects (via pylab or matplotlib). In particular, the pyplot.rcParams
    object can be used to scale the font size (rcParams["font.size"]) and line
    width (rcParams["lines.linewidth"]).
    :Parameters:
        label_func : callable
            A function to extract a label from a node. By default this is str(),
            but you can use a different function to select another string
            associated with each node. If this function returns None for a node,
            no label will be shown for that node.
        do_show : bool
            Whether to show() the plot automatically.
        show_confidence : bool
            Whether to display confidence values, if present on the tree.
        axes : matplotlib/pylab axes
            If a valid matplotlib.axes.Axes instance, the phylogram is plotted
            in that Axes. By default (None), a new figure is created.
        branch_labels : dict or callable
            A mapping of each clade to the label that will be shown along the
            branch leading to it. By default this is the confidence value(s) of
            the clade, taken from the ``confidence`` attribute, and can be
            easily toggled off with this function's ``show_confidence`` option.
            But if you would like to alter the formatting of confidence values,
            or label the branches with something other than confidence, then use
            this option.
        label_colors : dict or callable
            A function or a dictionary specifying the color of the tip label.
            If the tip label can't be found in the dict or label_colors is
            None, the label will be shown in black.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        try:
            import pylab as plt
        except ImportError:
            raise MissingPythonDependencyError(
                "Install matplotlib or pylab if you want to use draw."
            ) from None

    import matplotlib.collections as mpcollections

    # Arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []

    # Options for displaying branch labels / confidence
    def conf2str(conf):
        if int(conf) == conf:
            return str(int(conf))
        return str(conf)

    if not branch_labels:
        if show_confidence:

            def format_branch_label(clade):
                if clade.support is not None:
                    return conf2str(clade.support)
                return None
                # try:
                #     confidences = clade.supports
                #     # phyloXML supports multiple confidences
                # except AttributeError:
                #     pass
                # else:
                #     return "/".join(conf2str(cnf.value) for cnf in confidences)
                # if clade.confidence is not None:
                #     return conf2str(clade.confidence)
                # return None

        else:

            def format_branch_label(clade):
                return None

    elif isinstance(branch_labels, dict):

        def format_branch_label(clade):
            return branch_labels.get(clade)

    else:
        if not callable(branch_labels):
            raise TypeError(
                "branch_labels must be either a dict or a callable (function)"
            )
        format_branch_label = branch_labels

    # options for displaying label colors.
    if label_colors:
        if callable(label_colors):

            def get_label_color(label):
                return label_colors(label)

        else:
            # label_colors is presumed to be a dict
            def get_label_color(label):
                return label_colors.get(label, "black")

    else:

        def get_label_color(label):
            # if label_colors is not specified, use black
            return "black"

    # Layout

    def get_x_positions(tree):
        """Create a mapping of each clade to its horizontal position.
        Dict of {clade: x-coord}
        """
        # If there are no branch lengths, assume unit branch lengths
        depths = get_depths(tree)
        return depths

    def get_y_positions(tree):
        """Create a mapping of each clade to its vertical position.
        Dict of {clade: y-coord}.
        Coordinates are negative, and integers for tips.
        """
        maxheight = tree.count(tips=True)
        # Rows are defined by the tips
        heights = {
            tip: maxheight - i for i, tip in enumerate((tree.tips()))
        }

        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = (
                heights[clade.children[0]] + heights[clade.children[-1]]
            ) / 2.0

        if tree.is_root():
            calc_row(tree.root())
        return heights

    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)
    # The function draw_clade closes over the axes object
    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(1, 1, 1)
    elif not isinstance(axes, plt.matplotlib.axes.Axes):
        raise ValueError(f"Invalid argument for axes: {axes}")

    def draw_clade_lines(
        use_linecollection=False,
        orientation="horizontal",
        y_here=0,
        x_start=0,
        x_here=0,
        y_bot=0,
        y_top=0,
        color="black",
        lw=".1",
    ):
        """Create a line with or without a line collection object.
        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """
        if not use_linecollection and orientation == "horizontal":
            axes.hlines(y_here, x_start, x_here, color=color, lw=lw)
        elif use_linecollection and orientation == "horizontal":
            horizontal_linecollections.append(
                mpcollections.LineCollection(
                    [[(x_start, y_here), (x_here, y_here)]], color=color, lw=lw
                )
            )
        elif not use_linecollection and orientation == "vertical":
            axes.vlines(x_here, y_bot, y_top, color=color)
        elif use_linecollection and orientation == "vertical":
            vertical_linecollections.append(
                mpcollections.LineCollection(
                    [[(x_here, y_bot), (x_here, y_top)]], color=color, lw=lw
                )
            )

    def draw_clade(clade, x_start, color, lw):
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade]
        y_here = y_posns[clade]

        # Add node/taxon labels
        label = label_func(clade)
        label_color = get_label_color(label)
        if label not in (None, clade.__class__.__name__):
            # axes.text(
            #     x_here+1,
            #     y_here,
            #     f" {label}",
            #     verticalalignment="center",
            #     color=color,
            # )
            axes.annotate(
                f" {label} ",
                xy=(x_here,y_here),
                xytext=(-1,0),
                textcoords="offset points",
                horizontalalignment="right",
                verticalalignment="center",
                color=contrasting_color(label_color),
                backgroundcolor=label_color,
                bbox={'boxstyle':'Square, pad=0', 'color':label_color}
            )
            axes.plot(x_here, y_here, color=label_color, marker="s")


        # phyloXML-only graphics annotations
        if hasattr(clade, "color") and clade.color is not None:
            color = clade.color.to_hex()
        elif clade.is_tip():
            color = label_color
        else:
            color = 'k'

        if hasattr(clade, "width") and clade.width is not None:
            lw = clade.width * plt.rcParams["lines.linewidth"]
        # Draw a horizontal line from start to here
        draw_clade_lines(
            use_linecollection=True,
            orientation="horizontal",
            y_here=y_here,
            x_start=x_start,
            x_here=x_here,
            color=color,
            lw=lw,
        )
        if clade.children and len(clade.children) > 0:
            # Draw a vertical line connecting all children
            y_top = y_posns[clade.children[0]]
            y_bot = y_posns[clade.children[-1]]
            # Only apply widths to horizontal lines, like Archaeopteryx
            draw_clade_lines(
                use_linecollection=True,
                orientation="vertical",
                x_here=x_here,
                y_bot=y_bot,
                y_top=y_top,
                color=color,
                lw=lw,
            )
            # Draw descendents
            for child in clade:
                draw_clade(child, x_here, color, lw)

        # Add label above the branch (optional)
        conf_label = format_branch_label(clade)
        if conf_label:
            axes.text(
                0.5 * (x_start + x_here),
                y_here,
                conf_label,
                fontsize="small",
                horizontalalignment="center",
            )
    draw_clade(tree.root(), 0, "k", plt.rcParams["lines.linewidth"])

    # If line collections were used to create clade lines, here they are added
    # to the pyplot plot.
    for i in horizontal_linecollections:
        axes.add_collection(i)
    for i in vertical_linecollections:
        axes.add_collection(i)

    # Aesthetics

    try:
        name = tree.name
    except AttributeError:
        pass
    else:
        if name:
            axes.set_title(name)

    axes.get_xaxis().set_visible(False)
    axes.get_yaxis().set_visible(False)

    # Add margins around the tree to prevent overlapping the axes
    xmax = max(x_posns.values())
    axes.set_xlim(-0.05 * xmax, 1.25 * xmax)
    # Also invert the y-axis (origin at the top)
    # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
    axes.set_ylim(max(y_posns.values()) + 0.8, 0.2)

    axes.invert_xaxis()

    # Parse and process key word arguments as pyplot options
    for key, value in kwargs.items():
        try:
            # Check that the pyplot option input is iterable, as required
            list(value)
        except TypeError:
            raise ValueError(
                'Keyword argument "%s=%s" is not in the format '
                "pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict),"
                " or pyplot_option_name=(dict) " % (key, value)
            ) from None
        if isinstance(value, dict):
            getattr(plt, str(key))(**dict(value))
        elif not (isinstance(value[0], tuple)):
            getattr(plt, str(key))(*value)
        elif isinstance(value[0], tuple):
            getattr(plt, str(key))(*value[0], **dict(value[1]))

    if do_show:
        plt.show()




# 
# -----------------------------------------------------------------------------

def layout_tree(
    tree,
    label_func=str,
    show_confidence=True,
    # For power users
    branch_labels=None,
    clade_colors=None,
    base_clade_width=1,
    default_clade_color='#000000',
    default_label_color='#000000',
    identifier='name',
    *args,
    **kwargs,
):
    """Plot the given tree using matplotlib (or pylab).
    The graphic is a rooted tree, drawn with roughly the same algorithm as
    draw_ascii.
    Additional keyword arguments passed into this function are used as pyplot
    options. The input format should be in the form of:
    pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict), or
    pyplot_option_name=(dict).
    Example using the pyplot options 'axhspan' and 'axvline'::
        from Bio import Phylo, AlignIO
        from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
        constructor = DistanceTreeConstructor()
        aln = AlignIO.read(open('TreeConstruction/msa.phy'), 'phylip')
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        tree = constructor.upgma(dm)
        Phylo.draw(tree, axhspan=((0.25, 7.75), {'facecolor':'0.5'}),
        ... axvline={'x':0, 'ymin':0, 'ymax':1})
    Visual aspects of the plot can also be modified using pyplot's own functions
    and objects (via pylab or matplotlib). In particular, the pyplot.rcParams
    object can be used to scale the font size (rcParams["font.size"]) and line
    width (rcParams["lines.linewidth"]).
    :Parameters:
        label_func : callable
            A function to extract a label from a node. By default this is str(),
            but you can use a different function to select another string
            associated with each node. If this function returns None for a node,
            no label will be shown for that node.
        do_show : bool
            Whether to show() the plot automatically.
        show_confidence : bool
            Whether to display confidence values, if present on the tree.
        axes : matplotlib/pylab axes
            If a valid matplotlib.axes.Axes instance, the phylogram is plotted
            in that Axes. By default (None), a new figure is created.
        branch_labels : dict or callable
            A mapping of each clade to the label that will be shown along the
            branch leading to it. By default this is the confidence value(s) of
            the clade, taken from the ``confidence`` attribute, and can be
            easily toggled off with this function's ``show_confidence`` option.
            But if you would like to alter the formatting of confidence values,
            or label the branches with something other than confidence, then use
            this option.
        clade_colors : dict or callable
            A function or a dictionary specifying the color of the tip label.
            If the tip label can't be found in the dict or clade_colors is
            None, the label will be shown in black.
    """

    # Arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []
    nodes = []
    branch_nodes = []

    # Options for displaying branch labels / confidence
    def conf2str(conf):
        if int(conf) == conf:
            return str(int(conf))
        return str(conf)

    if not branch_labels:
        if show_confidence:

            def format_branch_label(clade):
                if clade.support is not None:
                    return conf2str(clade.support)
                return None

        else:

            def format_branch_label(clade):
                return None

    elif isinstance(branch_labels, dict):

        def format_branch_label(clade):
            return branch_labels.get(clade)

    else:
        if not callable(branch_labels):
            raise TypeError(
                "branch_labels must be either a dict or a callable (function)"
            )
        format_branch_label = branch_labels

    # options for displaying label colors.
    if clade_colors:
        if callable(clade_colors):

            def get_clade_color(label):
                return clade_colors(label)

        elif isinstance(clade_colors, dict):
            # clade_colors is presumed to be a dict
            def get_clade_color(label):
                return clade_colors.get(label, default_label_color)
        else:
            raise TypeError("clade_colors should be dict, callable, or False")

    else:

        def get_clade_color(label):
            # if clade_colors is not specified, use black
            return default_label_color

    # Layout

    def get_x_positions(tree):
        """Create a mapping of each clade to its horizontal position.
        Dict of {clade: x-coord}
        """
        # If there are no branch lengths, assume unit branch lengths
        depths = get_depths(tree)
        return depths

    def get_y_positions(tree):
        """Create a mapping of each clade to its vertical position.
        Dict of {clade: y-coord}.
        Coordinates are negative, and integers for tips.
        """
        maxheight = tree.count(tips=True)
        # Rows are defined by the tips
        heights = {
            tip: maxheight - i for i, tip in enumerate((tree.tips()))
        }

        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = (
                heights[clade.children[0]] + heights[clade.children[-1]]
            ) / 2.0

        if tree.is_root():
            calc_row(tree.root())
        return heights

    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)
    # The function draw_clade closes over the axes object

    def draw_clade_lines(
        orientation="horizontal",
        y_here=0,
        x_start=0,
        x_here=0,
        y_bot=0,
        y_top=0,
        color=default_clade_color,
        lw=".1",
    ):
        """Create a line with or without a line collection object.
        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """
        if orientation == "horizontal":
            horizontal_linecollections.append(
                dict(
                    x=x_start, y=y_here, x2=x_here, y2=y_here, color=color, lw=lw
                )
            )
        elif orientation == "vertical":
            vertical_linecollections.append(
                dict(
                    x=x_here, y=y_bot, x2=x_here, y2=y_top, color=color, lw=lw
                )
            )

    def draw_clade(clade, x_start, color, lw):
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade]
        y_here = y_posns[clade]

        # Add node/taxon labels
        label = label_func(clade)
        label_color = get_clade_color(clade)
        if label not in (None, clade.__class__.__name__):

            nodes.append(
                dict(
                    x=x_here,
                    y=y_here,
                    label=label,
                    color=(contrasting_color(label_color)
                           if label_color is not None else default_clade_color),
                    fill=label_color,
                    is_tip=clade.is_tip(),
                    **{identifier: clade.name}
                )
            )

        # phyloXML-only graphics annotations
        if hasattr(clade, "color") and clade.color is not None:
            color = clade.color.to_hex()
        elif clade.is_tip():
            color = label_color
        else:
            color = default_clade_color

        if hasattr(clade, "width") and clade.width is not None:
            lw = clade.width * base_clade_width

        # Draw a horizontal line from start to here
        draw_clade_lines(
            orientation="horizontal",
            y_here=y_here,
            x_start=x_start,
            x_here=x_here,
            color=color,
            lw=lw,
        )
        if clade.children and len(clade.children) > 0:
            # Draw a vertical line connecting all children
            y_top = y_posns[clade.children[0]]
            y_bot = y_posns[clade.children[-1]]
            # Only apply widths to horizontal lines, like Archaeopteryx
            draw_clade_lines(
                orientation="vertical",
                x_here=x_here,
                y_bot=y_bot,
                y_top=y_top,
                color=color,
                lw=lw,
            )
            # Draw descendents
            for child in clade:
                draw_clade(child, x_here, color, lw)

        # Add label above the branch (optional)
        conf_label = format_branch_label(clade)
        if conf_label:
            branch_nodes.append(
                dict(
                    x=0.5 * (x_start + x_here),
                    y=y_here,
                    label=conf_label
                )
            )
    draw_clade(tree.root(), 0, default_clade_color, base_clade_width)

    # If line collections were used to create clade lines, here they are added
    # to the pyplot plot.
    links = horizontal_linecollections + vertical_linecollections

    return dict(links=links, nodes=nodes)

    # Aesthetics

    # Add margins around the tree to prevent overlapping the axes
    # xmax = max(x_posns.values())
    # axes.set_xlim(-0.05 * xmax, 1.25 * xmax)
    # # Also invert the y-axis (origin at the top)
    # # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
    # axes.set_ylim(max(y_posns.values()) + 0.8, 0.2)

    # axes.invert_xaxis()


def fortify_tree(feature_data, tree, identifier=None, features=None, clade_colors=True, **kwargs):
    import pandas as pd
    from .utils import pretty_hex
    label_func = lambda c: pretty_hex(c.name) if (c.is_tip() and (c.name is not None)) else None


    # what should the column containing `tree_node.name` be called?
    # used both for joining to `feature_data` and for subsequent modifications
    # to the vega spec (e.g. tooltip, selection, etc.)
    if identifier is None:
        if feature_data.index.name is not None:
            identifier = feature_data.index.name
        else:
            identifier = 'feature'

    if features is not None:
        tip_names = [t.name for t in tree.tips()]
        tree = tree.shear(set(features).intersection(tip_names))

    if clade_colors == True:
        clade_colors = lambda c: hash_to_color(c.name) if c.name is not None else None
    
    # make layout giving x,y position of tree nodes and links
    layout = layout_tree(
        tree,
        label_func=label_func,
        clade_colors=clade_colors,
        show_confidence=False,
        identifier=identifier,
        **kwargs
    )

    nodes_df = pd.DataFrame(layout['nodes'])
    if feature_data is not None:
        nodes_df = nodes_df.join(feature_data, on=identifier)
        # if feature_data.index.name != identifier:
        #     feature_data = feature_data.set_index(identifier)
        # nodes_df = nodes_df.combine_first(feature_data)
        

    links_df = pd.DataFrame(layout['links'])

    return nodes_df, links_df


def plot_tree_alt(nodes_df, links_df, identifier='feature', label='hex', feature_scale=None, nodes={}, links={}):
    import altair as alt
    import pandas as pd
    
    if label == 'hex':
        from .utils import pretty_hex
        # label_func = lambda c: pretty_hex(c.name) if (c.is_tip() and (c.name is not None)) else None
        nodes_df['label'] = nodes_df[identifier].apply(
            lambda c: pretty_hex(c) if not pd.isna(c) else ''
        )
    elif label == 'mn':
        from .utils import hash_to_mn_short
        # label_func = lambda c: hash_to_mn_short(c.name) if (c.is_tip() and (c.name is not None)) else None
        nodes_df['label'] = nodes_df[identifier].apply(
            lambda c: hash_to_mn_short(c) if not pd.isna(c) else ''
        )
    elif callable(label):
        nodes_df['label'] = nodes_df[identifier].apply(label)


    if feature_scale is not None:
        _fill = alt.Fill('feature', scale=feature_scale)
    else:
        _fill = alt.Fill('color', scale=None)

    nodes_base = alt.Chart(nodes_df).encode(**{
        **dict(
            x=alt.X('x', axis=None),
            y=alt.Y('y', axis=None),
            # color=alt.Color('fill', scale=None),
            # fill=alt.Fill('fill',scale=None)
            fill=_fill,
            color=alt.value(None)
        ),
        **nodes
    })
    nodes_chart = nodes_base.mark_point()
    
    if label is not None:
        nodes_chart += nodes_base.mark_text(align='left', dx=5).encode(text='label')

    links_chart = alt.Chart(links_df).encode(**{
        **dict(
            x=alt.X('x',axis=None),
            x2='x2',
            y=alt.Y('y',axis=None),
            y2='y2', 
            # color=alt.Color('color',scale=None),
        ),
        **links
    }).mark_rule()

    chart = ((links_chart + nodes_chart)
        # .configure_legend(disable=True)
        # .configure_axis(grid=False, labels=False, title=None)
        # .configure_view(strokeOpacity=0)
    )
    return chart
