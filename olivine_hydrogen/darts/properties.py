"""Conservative H2/H2O/NaCl flash for an assumed isothermal packed-bed case.

H2 solubility: Chabab et al. (2020), doi:10.1016/j.ijhydene.2020.08.192.
The published dissolved mole fraction excludes salt from its denominator.
Water and NaCl are nonvolatile here; aqueous density/viscosity are declared
engineering inputs, not a brine EOS. No mineral equilibrium is implied.
"""
from pathlib import Path
import json
import numpy as np
from darts.physics.base.property_container import PropertyContainer
from darts.physics.properties.basic import ConstFunc, RockCompactionEvaluator

COMPONENTS=['H2','H2O','NaCl']
MW=np.array([2.01588,18.01528,58.4428])  # kg/kmol
FE_MW=55.845
G=9.80665e-5  # bar/m per kg/m3


def load_config(parameters=None):
    c=json.loads(Path(__file__).with_name('config.json').read_text())
    if isinstance(parameters,(str,Path)):parameters=json.loads(Path(parameters).read_text())
    if parameters:
        if set(parameters)-set(c):raise ValueError('Unknown config keys: '+str(set(parameters)-set(c)))
        c.update(parameters)
    for k in ['nx','ny','nz']:
        if not isinstance(c[k],int) or c[k]<2:raise ValueError('Grid dimensions must be integers >=2')
    for k,v in c.items():
        if isinstance(v,(int,float)) and not isinstance(v,bool) and (not np.isfinite(v) or v<=0):
            raise ValueError(k+' must be finite and positive')
    for k in ['fe_mass_fraction','ferrous_fraction','accessible_fraction','porosity','collection_efficiency','fuel_cell_efficiency']:
        if not 0<c[k]<=(1 if k!='porosity' else .99):raise ValueError('Invalid fraction '+k)
    if not 323.15<=c['temperature_k']<=373.15 or not 10<c['producer_bhp_bar']<230:
        raise ValueError('Solubility fit requires 323.15–373.15K and 10–230bar')
    if not 0<c['salinity_molal']<=5:raise ValueError('Salinity must lie in (0,5] mol/kg water')
    if not 10<=c['separator_pressure_bar']<=c['producer_bhp_bar']:raise ValueError('Separator pressure must lie between10bar and producerBHP')
    return c


def solubility_x(pressure_bar,temperature_k,molality):
    """Salt-free aqueous xH2, valid 10–230bar,323.15–373.15K,0–5molal.

    Tiny interpolation/Newton excursions are evaluated using the same smooth
    expression; the runner independently guards every accepted-state domain.
    No extra fugacity correction is applied to this empirical fit.
    """
    p=np.asarray(pressure_bar);t=temperature_k;m=np.asarray(molality)
    x0=p*(3.338844e-7*t+.0363161/t-.00020734)-2.1301815e-9*p*p
    return x0*np.exp(.018519*m*m-.30185103*m)


def hydrogen_density(p,t):
    """NIST pressure-explicit H2 density correlation, Lemmon et al.2008."""
    a=np.array([.05888460,-.06136111,-.002650473,.002731125,.001802374,-.001150707,.9588528e-4,-.1109040e-6,.1264403e-9])
    b=np.array([1.325,1.87,2.5,2.8,2.938,3.14,3.37,3.75,4.])
    d=np.array([1.,1.,2.,2.,2.42,2.63,3.,4.,5.])
    z=1+np.sum(a*(100/t)**b*(p/10)**d)
    return p*1e5*(MW[0]/1000)/(z*8.314462618*t)


class BrineFlash:
    def __init__(self,cfg):self.cfg=cfg
    def split(self,p,z):
        z=np.maximum(np.asarray(z,float),1e-14);z/=z.sum()
        h,w,s=z
        m=1000*s/(w*MW[1])
        xs=float(solubility_x(p,self.cfg['temperature_k'],m))
        if not 0<xs<1:raise ValueError('Invalid solubility outside property domain')
        h_diss=min(h,xs/(1-xs)*w)
        gas=max(h-h_diss,0.)
        nu=np.array([gas,1-gas]);x=np.zeros((2,3))
        x[0]=[1,0,0];x[1]=[h_diss,w,s];x[1]/=x[1].sum()
        if np.max(np.abs(nu@x-z))>1e-12:raise RuntimeError('Analytical flash did not conserve components')
        return nu,x
    def state(self,p,z):
        nu,x=self.split(p,z);t=self.cfg['temperature_k']
        rho=np.array([hydrogen_density(p,t),self.cfg['brine_density_kg_m3']*(1+self.cfg['brine_compressibility_per_bar']*(p-self.cfg['producer_bhp_bar']))])
        rm=rho/(x@MW);vol=nu/rm;sat=vol/vol.sum()
        return {'nu':nu,'x':x,'rho':rho,'rho_m':rm,'sat':sat,
                'molality':1000*x[1,2]/(x[1,1]*MW[1])}


class Density:
    def __init__(self,cfg,gas):self.cfg,self.gas=cfg,gas
    def evaluate(self,p,t,x):
        return hydrogen_density(p,t) if self.gas else self.cfg['brine_density_kg_m3']*(1+self.cfg['brine_compressibility_per_bar']*(p-self.cfg['producer_bhp_bar']))


class RelPerm:
    def evaluate(self,s):return float(np.clip(s,0,1)**2)


class Properties(PropertyContainer):
    def __init__(self,cfg):
        self.cfg=cfg;self.system=BrineFlash(cfg)
        super().__init__(['gas','brine'],COMPONENTS,MW,eps_z=1e-12,rock_comp=0.,temperature=cfg['temperature_k'])
        self.rock_compr_ev=RockCompactionEvaluator(pref=cfg['producer_bhp_bar'],compres=0.)
        self.density_ev={'gas':Density(cfg,True),'brine':Density(cfg,False)}
        self.viscosity_ev={'gas':ConstFunc(.0094),'brine':ConstFunc(cfg['brine_viscosity_cp'])}
        self.rel_perm_ev={'gas':RelPerm(),'brine':RelPerm()}
        self.capillary_pressure_ev=ConstFunc(np.zeros(2))
    def run_flash(self,pressure,temperature,zc,evaluate_PT=False):
        self.nu,self.x=self.system.split(pressure,zc)
        return np.flatnonzero(self.nu>0)
    def evaluate_mass_source(self,pressure,temperature,zc):
        # Reference reaction:3FeO+H2O -> Fe3O4+H2. Positive native KIN is
        # a sink; mesh.kin_factor supplies kmol extent/(m3 bulk day) relative
        # to this unit reference. All solids remain in the separate Fe ledger.
        self.mass_source[:]=[-1.,1.,0.]
        return self.mass_source
