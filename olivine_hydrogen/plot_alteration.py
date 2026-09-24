"""Conceptual mineral-volume sensitivity, not a coupled DARTS prediction."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch


def plot(output=None):
    output=Path(output or Path(__file__).with_name('results'))
    output.mkdir(parents=True,exist_ok=True)
    with plt.rc_context({'font.size':14,'axes.spines.top':False,'axes.spines.right':False}):
        fig,(ax,flow)=plt.subplots(1,2,figsize=(10,5.6),layout='constrained')
        for expansion,color in [(.2,'#397ca8'),(.4,'#609276'),(.6,'#b67642')]:
            limit=min(1.,.35/(.65*expansion))
            fraction=np.linspace(0,limit,160)
            porosity=.35-.65*expansion*fraction
            ax.plot(fraction,porosity,color=color,lw=2,label=f'Solid expansion {expansion:.0%}')
        ax.set(xlim=(0,1),ylim=(0,.37),xlabel='Altered mineral fraction, f',
               ylabel='Pore-volume fraction, φ',title='Fixed bulk volume:\npore filling only')
        ax.legend(fontsize=12,loc='upper right');ax.grid(alpha=.2)
        ax.text(.03,.025,'Initial porosity 0.35\nNo fracture, compaction or solid export',
                transform=ax.transAxes,fontsize=12)
        flow.set(xlim=(0,1),ylim=(0,1));flow.axis('off')
        flow.set_title('Competing effects\nduring alteration')
        boxes=[(.5,.90,'Olivine + water','#e6edde'),(.5,.67,'Hydrated secondary minerals\nand changed solid volume','#e7e9e9'),
               (.24,.38,'Stress and cracking\nif confinement and\nreaction allow it','#e1ebf1'),
               (.78,.38,'Pore filling and\nsurface coatings','#eee4da'),
               (.24,.10,'Fresh surfaces and\nnew fluid pathways','#e1ebf1'),
               (.78,.10,'Reduced fluid access\nand possible slowing','#eee4da')]
        for x,y,text,color in boxes:
            flow.text(x,y,text,ha='center',va='center',fontsize=12.5,
                      bbox={'boxstyle':'round,pad=0.45','facecolor':color,'edgecolor':'#778080'})
        for a,b in [((.5,.85),(.5,.74)),((.43,.62),(.24,.48)),((.59,.62),(.78,.47)),
                    ((.24,.28),(.24,.18)),((.78,.30),(.78,.18))]:
            flow.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=12,color='#555555'))
        fig.suptitle('Olivine alteration: geometric sensitivity and possible feedbacks',fontsize=16)
        fig.savefig(output/'alteration_feedbacks.png',dpi=190)
        plt.close(fig)
    return output/'alteration_feedbacks.png'


if __name__=='__main__':print(plot())
