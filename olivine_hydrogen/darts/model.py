"""Assumed 3-D reactive packed bed; this is not open-cavern Darcy flow."""
from __future__ import annotations
import numpy as np
from darts.engines import value_vector,well_control_iface
from darts.models.darts_model import DartsModel
from darts.nonlinear_solvers import NewtonSolver
from darts.physics.base.physics import PhysicsBase
from darts.reservoirs.struct_reservoir import StructReservoir
from .properties import COMPONENTS,MW,FE_MW,G,Properties,load_config


class Model(DartsModel):
    def __init__(self,parameters=None):
        self.cfg=load_config(parameters);c=self.cfg
        super().__init__();self.timer.node['initialization'].start()
        self.bulk_volume=c['rock_mass_kg']/c['rock_density_kg_m3']/(1-c['porosity'])
        self.height=self.bulk_volume/(c['length_x_m']*c['length_y_m'])
        self.nb=c['nx']*c['ny']*c['nz'];self.cell_volume=self.bulk_volume/self.nb
        self.depths=np.repeat((np.arange(c['nz'])+.5)*self.height/c['nz'],c['nx']*c['ny'])
        self.fe_initial_kmol=c['rock_mass_kg']*c['fe_mass_fraction']*c['ferrous_fraction']*c['accessible_fraction']/FE_MW
        self.fe_remaining=np.full(self.nb,self.fe_initial_kmol/self.nb)
        self.generated_kmol=0.;self._pending_extent=np.zeros(self.nb)
        self.reservoir=StructReservoir(self.timer,nx=c['nx'],ny=c['ny'],nz=c['nz'],
            dx=c['length_x_m']/c['nx'],dy=c['length_y_m']/c['ny'],dz=self.height/c['nz'],
            permx=c['permeability_md'],permy=c['permeability_md'],permz=c['permeability_md'],
            poro=c['porosity'],depth=self.depths,cache=False)
        self.property_container=Properties(c);self.system=self.property_container.system
        self.physics=PhysicsBase(COMPONENTS,['gas','brine'],self.timer,
            axes_step=[c['obl_pressure_step_bar'],c['obl_hydrogen_step'],c['obl_water_step']],
            axes_origin=[1.,1e-12,1e-12],epsilon_z=1e-12,extrapolation_flag=False,
            state_spec=PhysicsBase.StateSpecification.P,cache=False)
        self.physics.add_property_region(self.property_container)
        self.water_rate=c['water_feed_kg']/c['runtime_days']
        water_moles=self.water_rate/MW[1];salt_moles=self.water_rate*c['salinity_molal']/1000
        self.inj_composition=np.array([1e-10,(1-1e-10)*water_moles/(water_moles+salt_moles),(1-1e-10)*salt_moles/(water_moles+salt_moles)])
        self.inj_molar_rate=water_moles+salt_moles
        self.timer.node['initialization'].stop()
    def set_wells(self):
        c=self.cfg
        for name,ijk in [('INJ',(1,1,1)),('PROD',(c['nx'],c['ny'],c['nz']))]:
            self.reservoir.add_well(name)
            self.reservoir.add_perforation(name,res_cell_idx=ijk,well_diameter=2*c['well_radius_m'],well_indexD=0.)
    def set_well_controls(self):
        for w in self.reservoir.wells:
            if w.name=='INJ':
                self.physics.set_well_controls(wctrl=w.control,control_type=well_control_iface.MOLAR_RATE,is_inj=True,
                    target=self.inj_molar_rate,phase_name='brine',inj_composition=self.inj_composition.tolist())
            else:self.physics.set_well_controls(wctrl=w.control,control_type=well_control_iface.BHP,is_inj=False,target=self.cfg['producer_bhp_bar'])
    def set_initial_conditions(self):
        c=self.cfg;p=c['producer_bhp_bar']+G*c['brine_density_kg_m3']*(self.depths-self.depths[-1])
        water=1/(1+c['salinity_molal']*MW[1]/1000)
        self.physics.set_initial_conditions_from_array(self.reservoir.mesh,{'pressure':p,'H2':np.full(self.nb,1e-10),'H2O':np.full(self.nb,water*(1-1e-10))})
        self.reservoir.mesh.kin_factor=value_vector(np.zeros(self.nb))
    def set_solver(self):
        super().set_solver();c=self.cfg
        self.ts_control.dt_first=c['first_timestep_days'];self.ts_control.dt_min=1e-12
        self.ts_control.dt_max=c['max_timestep_days'];self.ts_control.dt_mult=1.5
        self.ts_control.runtime=c['runtime_days']
        self.nonlinear_solver=NewtonSolver(tolerance=c['nonlinear_tolerance'],max_iterations=20)
        self.linear_solver.spec.tolerance=1e-11
    def run_timestep(self,dt,t,verbose=0):
        """Freeze a conservative source over each attempted timestep.

        OBL source operators are immutable; the per-cell native kinetic factor
        carries the rate. Rejected attempts do not commit mineral consumption.
        """
        fraction=-np.expm1(-np.log(2)*dt/self.cfg['half_time_days']) if self.cfg['reaction_enabled'] else 0.
        extent=self.fe_remaining/3*fraction
        # Water limitation uses present cell inventory, never mineral capacity alone.
        q=self.cell_properties();water=self.cell_volume*self.cfg['porosity']*np.sum(q['sat']*q['rho_m']*q['x'][:,:,1],axis=1)
        self._pending_extent=np.minimum(extent,.99*water)
        self.reservoir.mesh.kin_factor=value_vector(self._pending_extent/(self.cell_volume*dt))
        return super().run_timestep(dt,t,verbose)
    def after_converged_timestep(self):
        self.fe_remaining-=3*self._pending_extent
        self.generated_kmol+=float(self._pending_extent.sum())
        if self.fe_remaining.min()<-1e-12:raise RuntimeError('Finite iron inventory became negative')
        super().after_converged_timestep()
    def reservoir_state(self):return np.asarray(self.physics.engine.X)[:3*self.nb].reshape(self.nb,3).copy()
    def cell_properties(self,states=None):
        states=self.reservoir_state() if states is None else states
        q=[self.system.state(float(s[0]),np.r_[s[1:],1-s[1:].sum()]) for s in states]
        return {k:np.array([v[k] for v in q]) for k in q[0]}
    def inventory(self,interpolated=False):
        states=self.reservoir_state();pv=self.cell_volume*self.cfg['porosity']
        if interpolated:
            ev=self.physics.reservoir_operators[0];kmol=np.zeros(3)
            for s in states:
                values=value_vector(np.zeros(ev.n_ops));self.physics.acc_flux_itor[0].evaluate(value_vector(s),values)
                kmol+=pv*np.asarray(values)[ev.ACC_OP:ev.ACC_OP+3]
            return {'component_kmol':kmol,'component_kg':kmol*MW}
        q=self.cell_properties(states);phase=pv*(q['sat']*q['rho_m'])[:,:,None]*q['x']
        return {'component_kmol':phase.sum(axis=(0,1)),'component_kg':phase.sum(axis=(0,1))*MW,
                'free_H2_kg':float(phase[:,0,0].sum()*MW[0]),'dissolved_H2_kg':float(phase[:,1,0].sum()*MW[0])}
    def perforation_rates(self):
        """Native TPFA flux; positive from reservoir to each well, kg/day.

        Uses interpolated FLUX×LAMBDA and actual native connection gravity/WI.
        Legacy engine.time_data component rates omit phase mobility here.
        """
        mesh=self.reservoir.mesh;states=np.asarray(self.physics.engine.X).reshape(mesh.n_blocks,3)
        bm,bp=np.asarray(mesh.block_m),np.asarray(mesh.block_p);trans,grav=np.asarray(mesh.tran),np.asarray(mesh.grav_coef)
        ev=self.physics.reservoir_operators[0];out={}
        for w in self.reservoir.wells:
            phase=np.zeros((2,3))
            for seg,cell,_wi,_wid in w.perforations:
                wc=w.well_body_idx+seg;idx=np.flatnonzero((bm==cell)&(bp==wc))
                if len(idx)!=1:raise RuntimeError('Unexpected well connection topology')
                k=idx[0];a=value_vector(np.zeros(ev.n_ops));b=value_vector(np.zeros(ev.n_ops))
                self.physics.acc_flux_itor[0].evaluate(value_vector(states[cell]),a)
                self.physics.acc_flux_w_itor.evaluate(value_vector(states[wc]),b)
                res,wel=np.asarray(a),np.asarray(b)
                for ph in range(2):
                    dp=states[wc,0]-states[cell,0]+.5*(res[ev.GRAV_OP+ph]+wel[ev.GRAV_OP+ph])*grav[k]-wel[ev.PC_OP+ph]+res[ev.PC_OP+ph]
                    up=res if dp<0 else wel
                    flow=-trans[k]*dp*up[ev.LAMBDA_OP+ph]
                    phase[ph]+=flow*up[ev.FLUX_OP+3*ph:ev.FLUX_OP+3*(ph+1)]
            out[w.name]={'component_kg_day':phase.sum(axis=0)*MW,'phase_component_kg_day':phase*MW}
        return out
