import numpy as np
from .ara import get_mcc
import dmritools as dm
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def label_rect(rect, text, ax, rotation=0):
    '''
    Utility function for labeling a Rectangle object.

    Parameters
    __________
    rect : matplotlib.patches.Rectangle
        Rectangle object to be labeled
    text : str
        Label text
    ax : matplotlib.axes._subplots.AxesSubplot
        Figure ax
    rotation : float
        Rotation angle for label
    '''
    ax.add_patch(rect)
    rx, ry = rect.get_xy()
    cx = rx + rect.get_width() / 2.0
    cy = ry + rect.get_height() / 2.0
    ax.annotate(text, (cx, cy), ha='center', va='center',
                rotation=rotation)


def connectome_viewer(data, cmap='inferno', ticks=[0, 1], ticklabels=[0, 1],
                      cmaptitle='Relative log\nweight'):
    '''
    Interface for generating connectome figure similar to "A mesoscale
    connectome of the mouse brain" (Oh et al, 2014).

    Parameters
    __________
    data : ndarray
        Connectome data
    cmap : str
        Matplotlib colormap
    ticks : list
        Colormap tick values
    ticklabels : list
        Colormap tick labels
    cmaptitle : str
        Colormap title

    Returns
    _______
    fig : matplotlib.figure.Figure
        Matplotlib figure
    ax : matplotlib.axes._subplots.AxesSubplot
        Figure ax
    '''

    # Assumes "data" does not include MDRN or fiber tracts

    allen_data = dm.get_connectome()
    mcc = get_mcc()
    tree = mcc.get_structure_tree()
    id_acro_map = tree.get_id_acronym_map()

    acros = ['Isocortex', 'OLF', 'HPF', 'CTXsp',
             'STR', 'PAL', 'TH', 'HY', 'MB', 'P', 'MY', 'CB']
    target_acros = [acro.split('-R')[0]
                    for acro in allen_data.columns[5:321] if acro.split('-')[0] not in ['MDRN', 'fiber tracts']]
    ny, nx = data.shape

    source_masks = {}
    target_masks = {}
    for acro in acros:
        source_masks[acro] = np.array([tree.structure_descends_from(
            id_acro_map[struct], id_acro_map[acro]) for struct in allen_data['primary-injection-structure']])
        target_masks[acro] = np.array([tree.structure_descends_from(
            id_acro_map[struct], id_acro_map[acro]) for struct in target_acros])

    colors = {'ipsi': '#d4d8de',
              'contra': '#c0c4cb',
              'source': '#c0c4cb',
              'Isocortex': '#97d3bb',
              'OLF': '#9dd08c',
              'HPF': '#07a56e',
              'CTXsp': '#b7dcb9',
              'STR': '#c7daf0',
              'PAL': '#7aa0d4',
              'TH': '#ee512c',
              'HY': '#f58e8d',
              'MB': '#8568ae',
              'P': '#f47836',
              'MY': '#f7e2e5',
              'CB':  '#f5ee74'}

    text = {'ipsi': 'Target: right hemisphere (ipsilateral)',
            'contra': 'Target: left hemisphere (contralateral)',
            'source': 'Source (injections)',
            'Isocortex': 'Iso-\ncortex',
            'OLF': 'OLF',
            'HPF': 'HPF',
            'CTXsp': 'CTXsp',
            'STR': 'STR',
            'PAL': 'PAL',
            'TH': 'Thal',
            'HY': 'Hypo-\nthal',
            'MB': 'Mid-\nbrain',
            'P': 'Pons',
            'MY': 'Medulla',
            'CB':  'CB'}

    hemilabel_dx = 20
    structlabel_dx = 40
    buff = hemilabel_dx + structlabel_dx
    mat = np.zeros((ny + buff, nx + buff))
    mat[buff:, buff:] = data

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_axis_off()

    im = ax.imshow(mat, cmap=cmap, vmin=ticks[0], vmax=ticks[1])

    ipsi_rect = patches.Rectangle(xy=(buff-0.5, -0.5), width=nx // 2,
                                  height=hemilabel_dx, facecolor=colors['ipsi'])
    contra_rect = patches.Rectangle(xy=(buff + nx//2 - 0.5, -0.5), width=nx // 2,
                                    height=hemilabel_dx, facecolor=colors['contra'])
    source_rect = patches.Rectangle(xy=(-0.5, buff-0.5), width=hemilabel_dx,
                                    height=ny, facecolor=colors['source'])
    label_rect(ipsi_rect, text['ipsi'], ax)
    label_rect(contra_rect, text['contra'], ax)
    label_rect(source_rect, text['source'], ax, rotation=90)

    corner_rect = patches.Rectangle(xy=(-0.5, -0.5), width=buff,
                                    height=buff, facecolor='white')
    ax.add_patch(corner_rect)

    for acro in acros:
        src_rect = ax.add_patch(patches.Rectangle(xy=(hemilabel_dx - 0.5, buff-0.5 + np.where(source_masks[acro])[0][0]),
                                                  width=structlabel_dx,
                                                  height=source_masks[acro].sum(
        ),
            facecolor=colors[acro]))
        label_rect(src_rect, text[acro], ax)

        if acro in ['OLF', 'HPF', 'CTXsp', 'STR', 'PAL', 'P']:
            rotation = 90
        else:
            rotation = 0

        tgt_r_rect = ax.add_patch(patches.Rectangle(xy=(buff - 0.5 + np.where(target_masks[acro])[0][0], hemilabel_dx-0.5),
                                                    width=target_masks[acro].sum(
        ),
            height=structlabel_dx,
            facecolor=colors[acro]))
        label_rect(tgt_r_rect, text[acro], ax, rotation=rotation)

        tgt_l_rect = ax.add_patch(patches.Rectangle(xy=(nx // 2 + buff - 0.5 + np.where(target_masks[acro])[0][0], hemilabel_dx-0.5),
                                                    width=target_masks[acro].sum(
        ),
            height=structlabel_dx,
            facecolor=colors[acro]))

        label_rect(tgt_l_rect, text[acro], ax, rotation=rotation)

    fig.tight_layout()
    cbax = inset_axes(ax, width="75%", height="20%", loc='center',
                      bbox_to_anchor=corner_rect.get_bbox(),
                      bbox_transform=ax.transData,
                      borderpad=0)
    cbax.set_title(cmaptitle, fontsize=8)

    cbar = fig.colorbar(im, cax=cbax, orientation='horizontal',
                        ticks=ticks)
    cbar.ax.set_xticklabels(ticklabels, fontdict={'fontsize': 6.5})
    cbar.ax.set_xlabel(r'log$_{10}$', fontsize=6.5, labelpad=-10)
    return fig, ax
