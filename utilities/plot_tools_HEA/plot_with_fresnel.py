#!/usr/bin/env python3

#ASIF
# https://docs.matador.science/en/latest/notebooks/interactive/visualisation.html

from matador.utils.viz_utils import fresnel_plot, fresnel_view
from matador.utils.optimade_utils import optimade2dict
import requests
import matplotlib.pyplot as plt
import ase.build
from ase.build import bulk
import fresnel
from functools import partial
from matador.crystal.network import draw_network
from matador.utils.ase_utils import ase2dict
from matador.crystal import Crystal
from matador.utils.viz_utils import rerender_scenes_to_axes
import PIL
from PIL import ImageShow

def plot(crystals):
    fig, axes, scenes = fresnel_plot(
        crystals,
        figsize=(15,15),
        fig_cols=1,
        fig_rows=1,
        extension=[(2, 3, 3)],
        lights=[fresnel.light.butterfly],
        standardize=False,
        renderer=partial(fresnel.pathtrace, w=1920, h=1080, light_samples=64, samples=64),
        #renderer=partial(fresnel.preview, w=1920, h=1080, anti_alias=True),
    )
    tt = scenes
    plt.savefig('test.png', bbox_inches = 'tight', dpi= 300)


cu = bulk('Fe', 'bcc', a=2.832, orthorhombic=True)
t = Crystal.from_ase(cu)

print( (t))

plot(t)

