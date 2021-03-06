'''
Script with helper methods for style settings
'''

from ROOT import gStyle, TGaxis  # pylint: disable=import-error,no-name-in-module
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def formataxes(axes, xmin, xmax, ymin, ymax, xtitle, ytitle,
               majticksx=None, minticksx=None, majticksy=None, minticksy=None):
    axes.set_title('')
    axes.set_xlabel(xtitle, fontsize=14)
    axes.set_ylabel(ytitle, fontsize=14)
    axes.set_xlim(xmin, xmax)
    axes.set_ylim(ymin, ymax)
    majorFormatter = FormatStrFormatter('%.1f')
    if majticksx and minticksx:
        majorLocator_x = MultipleLocator(majticksx)
        minorLocator_x = MultipleLocator(minticksx)
        axes.xaxis.set_major_locator(majorLocator_x)
        axes.xaxis.set_minor_locator(minorLocator_x)
        axes.xaxis.set_major_formatter(majorFormatter)
    if majticksy and minticksy:
        majorLocator_y = MultipleLocator(majticksy)
        minorLocator_y = MultipleLocator(minticksy)
        axes.yaxis.set_major_locator(majorLocator_y)
        axes.yaxis.set_minor_locator(minorLocator_y)
        axes.yaxis.set_major_formatter(majorFormatter)
    axes.xaxis.set_ticks_position('both')
    axes.yaxis.set_ticks_position('both')
    axes.tick_params(length=8, width=1)
    axes.tick_params(length=4, width=1, which='minor')
    axes.set_rasterization_zorder(1)


# pylint: disable=too-many-branches
def SetGlobalStyle(**kwargs):
    '''
    Method to set global style.

    Parameters
    ----------

    - padrightmargin (float), default = 0.035
    - padleftmargin (float), default = 0.12
    - padtopmargin (float), default = 0.035
    - padbottommargin (float), default = 0.12

    - titlesize (float), default = 0.050
    - titlesizex (float), default = 0.050
    - titlesizey (float), default = 0.050
    - titlesizez (float), default = 0.050

    - labelsize (float), default = 0.045
    - labelsizex (float), default = 0.045
    - labelsizey (float), default = 0.045
    - labelsizez (float), default = 0.045

    - titleoffset (float), default = 1.2
    - titleoffsetx (float), default = 1.2
    - titleoffsey (float), default = 1.2
    - titleoffsetz (float), default = 1.2

    - opttitle (int), default = 0
    - optstat (int), default = 0

    - maxdigits (int), default no max value
    '''

    # pad margins
    if 'padrightmargin' in kwargs:
        gStyle.SetPadRightMargin(kwargs['padrightmargin'])
    else:
        gStyle.SetPadRightMargin(0.035)

    if 'padleftmargin' in kwargs:
        gStyle.SetPadLeftMargin(kwargs['padleftmargin'])
    else:
        gStyle.SetPadLeftMargin(0.12)

    if 'padtopmargin' in kwargs:
        gStyle.SetPadTopMargin(kwargs['padtopmargin'])
    else:
        gStyle.SetPadTopMargin(0.035)

    if 'padbottommargin' in kwargs:
        gStyle.SetPadBottomMargin(kwargs['padbottommargin'])
    else:
        gStyle.SetPadBottomMargin(0.1)

    # title sizes
    if 'titlesize' in kwargs:
        gStyle.SetTitleSize(kwargs['titlesize'], 'xyz')
    else:
        gStyle.SetTitleSize(0.050, 'xyz')

    if 'titlesizex' in kwargs:
        gStyle.SetTitleSize(kwargs['titlesizex'], 'x')
    if 'titlesizey' in kwargs:
        gStyle.SetTitleSize(kwargs['titlesizex'], 'y')
    if 'titlesizez' in kwargs:
        gStyle.SetTitleSize(kwargs['titlesizex'], 'z')

    # label sizes
    if 'labelsize' in kwargs:
        gStyle.SetLabelSize(kwargs['labelsize'], 'xyz')
    else:
        gStyle.SetLabelSize(0.045, 'xyz')

    if 'labelsizex' in kwargs:
        gStyle.SetLabelSize(kwargs['labelsizex'], 'x')
    if 'labelsizey' in kwargs:
        gStyle.SetLabelSize(kwargs['labelsizex'], 'y')
    if 'labelsizez' in kwargs:
        gStyle.SetLabelSize(kwargs['labelsizex'], 'z')

    # title offsets
    if 'titleoffset' in kwargs:
        gStyle.SetTitleOffset(kwargs['titleoffset'], 'xyz')
    else:
        gStyle.SetTitleOffset(1.2, 'xyz')

    if 'titleoffsetx' in kwargs:
        gStyle.SetTitleOffset(kwargs['titleoffsetx'], 'x')
    if 'titleoffsety' in kwargs:
        gStyle.SetTitleOffset(kwargs['titleoffsety'], 'y')
    if 'titleoffsetz' in kwargs:
        gStyle.SetTitleOffset(kwargs['titleoffsetz'], 'z')

    # other options
    if 'opttitle' in kwargs:
        gStyle.SetOptTitle(kwargs['opttitle'])
    else:
        gStyle.SetOptTitle(0)

    if 'optstat' in kwargs:
        gStyle.SetOptStat(kwargs['optstat'])
    else:
        gStyle.SetOptStat(0)

    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetLegendBorderSize(0)

    if 'maxdigits' in kwargs:
        TGaxis.SetMaxDigits(kwargs['maxdigits'])


def SetObjectStyle(obj, **kwargs):
    '''
    Method to set root object style.

    Parameters
    ----------

    - obj: object to set style

    - linecolor (int) default 1 (black)
    - linealpha (float) default 1
    - linewitdh (int) default 2
    - linestyle (int) default 1

    - markercolor (int) default 1 (black)
    - markeralpha (float) default 1
    - markerstyle (int) default 20 (full circle)
    - markersize (int) default 20 (full circle)

    - fillcolor (int) default no filling
    - fillalpha (float) default 1
    - fillstyle (int) default 0 (no style)

    - color (int) sets same color for line, marker and fill
    - alpha (float) sets same alpha for line, marker and fill
    '''

    # alpha parameters
    lalpha = kwargs.get('linealpha', 1)
    malpha = kwargs.get('markeralpha', 1)
    falpha = kwargs.get('fillalpha', 1)
    if 'alpha' in kwargs:
        lalpha = kwargs['linealpha']
        malpha = kwargs['markeralpha']
        falpha = kwargs['fillalpha']

    # line styles
    if 'linecolor' in kwargs:
        obj.SetLineColorAlpha(kwargs['linecolor'], lalpha)
    else:
        obj.SetLineColorAlpha(1, lalpha)

    if 'linewidth' in kwargs:
        obj.SetLineWidth(kwargs['linewidth'])
    else:
        obj.SetLineWidth(2)

    if 'linestyle' in kwargs:
        obj.SetLineStyle(kwargs['linestyle'])
    else:
        obj.SetLineStyle(1)

    # marker styles
    if 'markercolor' in kwargs:
        obj.SetMarkerColorAlpha(kwargs['markercolor'], malpha)
    else:
        obj.SetMarkerColorAlpha(1, malpha)

    if 'markersize' in kwargs:
        obj.SetMarkerSize(kwargs['markersize'])
    else:
        obj.SetMarkerSize(1)

    if 'markerstyle' in kwargs:
        obj.SetMarkerStyle(kwargs['markerstyle'])
    else:
        obj.SetMarkerStyle(20)

    # fill styles
    if 'fillcolor' in kwargs:
        obj.SetFillColorAlpha(kwargs['fillcolor'], falpha)

    if 'fillstyle' in kwargs:
        obj.SetFillStyle(kwargs['fillstyle'])

    #global color
    if 'color' in kwargs:
        obj.SetLineColorAlpha(kwargs['color'], lalpha)
        obj.SetMarkerColorAlpha(kwargs['color'], malpha)
        obj.SetFillColorAlpha(kwargs['color'], falpha)
