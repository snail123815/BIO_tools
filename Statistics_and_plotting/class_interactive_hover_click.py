class interactionHoverClick:
    '''
    Working with scatter plot, put label as mouse hovering and clicked
    Usage:

    HoverClick = interactionHoverClick(colorDict, fig, ax)

    fig.canvas.mpl_connect("motion_notify_event", HoverClick.hover)
    fig.canvas.mpl_connect("button_press_event", HoverClick.click)

    colorDict -> {ax0:pd.Series(index=data2plot.index, dtype=np.dtype('U')),
                  ax1:pd.Series(index=data2plot.index, dtype=np.dtype('U'))
                  }
    The indexing method is done through the numbering of the index. Thus, make sure
    the index is the same as the actuall ploted dataframe.

    '''

    def __init__(self, colorDict, fig, axs):
        self.colorDict = colorDict
        self.fig = fig
        self.axs = axs
        self.hovered = []

        # set one annotation (empty) for each plot then hide them
        self.annotAxs = self.creatHoverAnnotation(axs)

        self.clicked = self.containerClicked(axs)
    # __init__

    def containerClicked(self, axs):
        '''Make a dictionary storing all clicked dots, for
        later query and clear (if found) the annotation.
        '''
        clicked = {}
        try:
            for ax in self.axs.ravel():  # if multiple plots passed in.
                clicked[ax] = {}
        except:  # when axs.ravel() do not work, that means there is one plot
            ax = self.axs
            clicked[ax] = {}
        return clicked
    # containerClicked

    def creatHoverAnnotation(self, axs):
        '''Annotation style for mouse hover, will set one annotation (empty) for each plot
        then hide them.'''
        annotAxs = {}
        try:
            for ax in axs.ravel():  # same reason as fuc. containerClicked()
                annotAxs[ax] = ax.annotate('',
                                           (0, 0),
                                           xytext=(-40, -30),
                                           textcoords="offset points",
                                           bbox=dict(boxstyle='round', fc='w', alpha=0.8),
                                           arrowprops=dict(arrowstyle='->')
                                           )
                annotAxs[ax].set_visible(False)
        except:
            ax = axs
            annotAxs[ax] = ax.annotate('',
                                       (0, 0),
                                       xytext=(-40, -30),
                                       textcoords="offset points",
                                       bbox=dict(boxstyle='round', fc='w', alpha=0.8),
                                       arrowprops=dict(arrowstyle='->')
                                       )
            annotAxs[ax].set_visible(False)
        return annotAxs
    # creatHoverAnnotation

    def hover(self, event):
        if event.inaxes != None:
            ax = event.inaxes
            # class matplotlib.collections.PathCollection, here refere to the scatter plot
            sc = ax.get_children()[0]
            cont, ind = sc.contains(event)
            annot = self.annotAxs[ax]
            if cont:
                try:
                    exist = (ind['ind'][0] in self.hovered)
                except:
                    exist = False
                if not exist:
                    hovered = ind['ind'][0]
                    pos = sc.get_offsets()[ind['ind'][0]]
                    text = '\n'.join(self.colorDict[ax].index[i] for i in ind['ind'])
                    annot.xy = pos
                    annot.set_text(text)
                    annot.set_visible(True)
                    self.fig.canvas.draw_idle()
            else:
                if annot.get_visible():
                    annot.set_visible(False)
                    self.fig.canvas.draw_idle()
    # hover

    def click(self, event):
        if event.inaxes != None:
            ax = event.inaxes
            # class matplotlib.collections.PathCollection, here refere to the scatter plot
            sc = ax.get_children()[0]
            cont, ind = sc.contains(event)
            if cont:
                if ind['ind'][0] in self.clicked[ax]:
                    self.clicked[ax][ind['ind'][0]].remove()
                    self.clicked[ax].pop(ind['ind'][0])
                    self.fig.canvas.draw_idle()
                else:
                    self.clicked[ax][ind['ind'][0]] = ax.annotate('\n'.join(self.colorDict[ax].index[i] for i in ind['ind']),
                                                                  sc.get_offsets()[ind['ind'][0]],
                                                                  xytext=(-30, -20),
                                                                  textcoords="offset points",
                                                                  bbox=dict(boxstyle='round',
                                                                            fc='r', alpha=0.3),
                                                                  arrowprops=dict(arrowstyle='->'),
                                                                  size='xx-small'
                                                                  )

                    self.fig.canvas.draw_idle()
    # click
