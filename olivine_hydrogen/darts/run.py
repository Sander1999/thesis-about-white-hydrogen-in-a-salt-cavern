"""Run the assumed 3-D packed-bed alternative with native DARTS.

From project root: ./run_darts.sh olivine_hydrogen/darts/run.py --suite
"""
from __future__ import annotations
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
if str(Path(__file__).resolve().parents[2]) not in sys.path:sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from darts.engines import redirect_darts_output,set_num_threads,value_vector
from olivine_hydrogen.darts.model import Model
from olivine_hydrogen.darts.properties import COMPONENTS,MW,FE_MW


class RecordedModel(Model):
    def __init__(self,parameters=None):super().__init__(parameters);self.records=[]
    def record(self):
        t=float(self.physics.engine.t);state=self.reservoir_state();q=self.cell_properties()
        if not np.isfinite(state).all() or state[:,0].min()<10 or state[:,0].max()>230:
            raise RuntimeError('Accepted state outside 10–230bar solubility domain')
        if q['molality'].min()<0 or q['molality'].max()>5:raise RuntimeError('Accepted salinity outside 0–5molal')
        exact=self.inventory();native=self.inventory(True)
        rates=self.perforation_rates() if t else {n:{'component_kg_day':np.zeros(3),'phase_component_kg_day':np.zeros((2,3))} for n in ['INJ','PROD']}
        prev=self.records[-1] if self.records else None;dt=t-prev['time_days'] if prev else 0
        row={'time_days':t,'pressure_min_bar':float(state[:,0].min()),'pressure_max_bar':float(state[:,0].max()),
             'generated_H2_kg':self.generated_kmol*MW[0],'unreacted_accessible_Fe_kg':float(self.fe_remaining.sum()*FE_MW),
             'free_H2_kg':exact['free_H2_kg'],'dissolved_H2_kg':exact['dissolved_H2_kg'],
             'gas_H2_out_kg':(prev['gas_H2_out_kg'] if prev else 0)+dt*rates['PROD']['phase_component_kg_day'][0,0],
             'aqueous_H2_out_kg':(prev['aqueous_H2_out_kg'] if prev else 0)+dt*rates['PROD']['phase_component_kg_day'][1,0],
             'max_gas_saturation':float(q['sat'][:,0].max()),'salinity_min_molal':float(q['molality'].min()),'salinity_max_molal':float(q['molality'].max())}
        # Equilibrium flash of each instantaneous total exported stream at
        # separator pressure; do not mistake dissolved export for collected H2.
        molar_out=np.maximum(rates['PROD']['component_kg_day']/MW,0.)
        gas_rate=0.;outlet_ratio=0.
        if molar_out.sum()>0:
            nu,x=self.system.split(self.cfg['separator_pressure_bar'],molar_out/molar_out.sum())
            gas_rate=float(nu[0]*molar_out.sum()*MW[0])
            outlet_ratio=float(molar_out[0]/max(molar_out[1],1e-30))
        row['outlet_H2_mol_per_mol_water']=outlet_ratio
        row['separator_gas_H2_rate_kg_day']=gas_rate
        row['separator_H2_available_kg']=(prev['separator_H2_available_kg'] if prev else 0.)+dt*gas_rate
        row['collected_H2_kg']=row['separator_H2_available_kg']*self.cfg['collection_efficiency']
        row['gross_electricity_kwh']=row['collected_H2_kg']*self.cfg['fuel_cell_efficiency']*self.cfg['hydrogen_lhv_kwh_kg']
        reaction=np.array([1.,-1.,0.])*self.generated_kmol*MW
        for j,name in enumerate(COMPONENTS):
            row[name+'_kg']=float(exact['component_kg'][j]);row['native_'+name+'_kg']=float(native['component_kg'][j])
            row[name+'_in_kg']=(prev[name+'_in_kg'] if prev else 0)-dt*rates['INJ']['component_kg_day'][j]
            row[name+'_out_kg']=(prev[name+'_out_kg'] if prev else 0)+dt*rates['PROD']['component_kg_day'][j]
            first=self.records[0] if self.records else row
            ledger=row[name+'_out_kg']-row[name+'_in_kg']-reaction[j]
            row[name+'_balance_kg']=row['native_'+name+'_kg']-first['native_'+name+'_kg']+ledger
            row[name+'_exact_balance_kg']=row[name+'_kg']-first[name+'_kg']+ledger
        self.records.append(row)
    def after_converged_timestep(self):super().after_converged_timestep();self.record()


def write_csv(path,rows):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def coordinates(m):
    c=m.cfg
    z,y,x=np.meshgrid((np.arange(c['nz'])+.5)*m.height/c['nz'],(np.arange(c['ny'])+.5)*c['length_y_m']/c['ny'],(np.arange(c['nx'])+.5)*c['length_x_m']/c['nx'],indexing='ij')
    return np.column_stack((x.ravel(),y.ravel(),z.ravel()))


def internal_velocity(m):
    """Average native Darcy face velocity per axis, m/day; brine phase only."""
    mesh=m.reservoir.mesh;s=m.reservoir_state();xyz=coordinates(m);vel=np.zeros((m.nb,3));counts=np.zeros((m.nb,3))
    ev=m.physics.reservoir_operators[0];ops=[]
    for state in s:
        a=value_vector(np.zeros(ev.n_ops));m.physics.acc_flux_itor[0].evaluate(value_vector(state),a);ops.append(np.array(a))
    dims=np.array([m.cfg['length_x_m']/m.cfg['nx'],m.cfg['length_y_m']/m.cfg['ny'],m.height/m.cfg['nz']])
    for a,b,trans,grav in zip(mesh.block_m,mesh.block_p,mesh.tran,mesh.grav_coef):
        if a>=m.nb or b>=m.nb:continue
        ph=1;dp=s[b,0]-s[a,0]+.5*(ops[a][ev.GRAV_OP+ph]+ops[b][ev.GRAV_OP+ph])*grav
        up=ops[a] if dp<0 else ops[b];flow=-trans*dp*up[ev.LAMBDA_OP+ph]
        displacement=xyz[b]-xyz[a];axis=int(np.argmax(np.abs(displacement)))
        v=flow/(np.prod(dims)/dims[axis])*np.sign(displacement[axis])
        for cell in [a,b]:vel[cell,axis]+=v;counts[cell,axis]+=1
    return np.divide(vel,counts,out=np.zeros_like(vel),where=counts>0)


def run_case(parameters=None,output=None,*,plots=True):
    output=Path(output or Path(__file__).with_name('output')/'reactive').resolve();output.mkdir(parents=True,exist_ok=True)
    (output/'status.json').write_text('{"status":"running"}\n')
    try:
        redirect_darts_output(str(output/'native.log'));set_num_threads(2)
        m=RecordedModel(parameters);m.init();m.set_output(output_folder=str(output));m.record()
        m.run(m.cfg['runtime_days'],verbose=0)
        if abs(m.physics.engine.t-m.cfg['runtime_days'])>1e-10:raise RuntimeError('DARTS stopped early')
        rows=m.records;first,last=rows[0],rows[-1];balance={};exact_balance={}
        for name in COMPONENTS:
            scale=max(first[name+'_kg']+last[name+'_in_kg']+(last['generated_H2_kg'] if name=='H2' else 0),1e-9)
            balance[name]=max(abs(r[name+'_balance_kg']) for r in rows)/scale
            exact_balance[name]=max(abs(r[name+'_exact_balance_kg']) for r in rows)/scale
        capacity=m.fe_initial_kmol/3*MW[0]
        expected=capacity*(-np.expm1(-np.log(2)*m.cfg['runtime_days']/m.cfg['half_time_days'])) if m.cfg['reaction_enabled'] else 0
        summary={'scope':'Assumed 3D packed-bed flow-through alternative; uncalibrated reaction time and permeability. Not open-cavern Darcy flow.',
                 'parameters':m.cfg,'days':last['time_days'],'cells':m.nb,'height_m':m.height,'bulk_volume_m3':m.bulk_volume,
                 'pore_volume_m3':m.bulk_volume*m.cfg['porosity'],'initial_pore_water_kg':first['H2O_kg'],'initial_pore_salt_kg':first['NaCl_kg'],
                 'target_cumulative_water_feed_kg':m.cfg['water_feed_kg'],'actual_cumulative_water_in_kg':last['H2O_in_kg'],
                 'actual_cumulative_salt_in_kg':last['NaCl_in_kg'],'redox_capacity_H2_kg':capacity,'generated_H2_kg':last['generated_H2_kg'],
                 'remaining_accessible_Fe_kg':last['unreacted_accessible_Fe_kg'],'aqueous_H2_exported_kg':last['aqueous_H2_out_kg'],
                 'gas_H2_exported_kg':last['gas_H2_out_kg'],'separator_H2_available_kg':last['separator_H2_available_kg'],
                 'collected_H2_kg':last['collected_H2_kg'],'gross_electricity_kwh':last['gross_electricity_kwh'],
                 'max_outlet_H2_mol_per_mol_water':max(r['outlet_H2_mol_per_mol_water'] for r in rows),
                 'remaining_free_H2_kg':last['free_H2_kg'],'remaining_dissolved_H2_kg':last['dissolved_H2_kg'],
                 'reaction_water_consumed_kg':m.generated_kmol*MW[1],'solid_oxygen_gain_kg':m.generated_kmol*(MW[1]-MW[0]),
                 'analytical_generated_H2_kg':expected,'reaction_integration_error_kg':last['generated_H2_kg']-expected,
                 'max_native_relative_balance_error':balance,'max_exact_relative_balance_error':exact_balance,
                 'final_pressure_range_bar':[last['pressure_min_bar'],last['pressure_max_bar']],
                 'accepted_timesteps':len(rows)-1,'numerical_H2_initial_kg':first['H2_kg'],'numerical_H2_injected_kg':last['H2_in_kg'],
                 'balance_method':'Native OBL storage and independently reconstructed TPFA well fluxes, right-endpoint integration over accepted backward-Euler timesteps; solid source committed only on acceptance.',
                 'collection_notice':'Each instantaneous exported H2/water/salt stream is equilibrated at separator pressure and unchanged temperature, then collection and fuel-cell efficiencies are applied. This flow-through separator differs from the retained batch model. Gross electricity excludes pumping/heating/processing demand.'}
        if max(balance.values())>2e-5:raise RuntimeError('Native mass balance failed: '+str(balance))
        if abs(summary['reaction_integration_error_kg'])>max(capacity*1e-10,1e-10):raise RuntimeError('Finite-reaction integration failed')
        write_csv(output/'history.csv',rows);xyz=coordinates(m);q=m.cell_properties();velocity=internal_velocity(m)
        concentration=np.sum(q['sat']*q['rho_m']*q['x'][:,:,0],axis=1)*MW[0]
        fields=[{'cell':i,'x_m':xyz[i,0],'y_m':xyz[i,1],'depth_m':xyz[i,2],'pressure_bar':m.reservoir_state()[i,0],
                 'H2_kg_m3_pore_fluid':concentration[i],'gas_saturation':q['sat'][i,0],'aqueous_H2_salt_free_mole_fraction':q['x'][i,1,0]/(q['x'][i,1,0]+q['x'][i,1,1]),
                 'salinity_molal':q['molality'][i],'darcy_x_m_day':velocity[i,0],'darcy_y_m_day':velocity[i,1],'darcy_z_m_day':velocity[i,2]} for i in range(m.nb)]
        write_csv(output/'spatial_final.csv',fields)
        (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        if plots:plot_case(m,rows,xyz,concentration,velocity,output)
        (output/'status.json').write_text('{"status":"complete"}\n')
        return summary,rows
    except Exception as error:
        (output/'status.json').write_text(json.dumps({'status':'failed','error':str(error)})+'\n');raise


def plot_case(m,rows,xyz,concentration,velocity,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    t=np.array([r['time_days'] for r in rows])*24
    fig,axes=plt.subplots(1,3,figsize=(13,4),layout='constrained')
    for key,label in [('generated_H2_kg','Generated'),('aqueous_H2_out_kg','Exported in brine'),('dissolved_H2_kg','Dissolved in bed'),('gas_H2_out_kg','Exported gas'),('collected_H2_kg','Collected after separator')]:
        axes[0].plot(t,[r[key]*1000 for r in rows],label=label)
    axes[0].set_ylabel('Hydrogen (g)');axes[0].legend(fontsize=8)
    axes[1].plot(t,[r['H2O_in_kg']/1000 for r in rows],label='Water feed');axes[1].plot(t,[r['NaCl_in_kg']/1000 for r in rows],label='Salt feed');axes[1].set_ylabel('Cumulative feed (tonne)');axes[1].legend()
    axes[2].plot(t,[r['pressure_min_bar'] for r in rows],label='Minimum');axes[2].plot(t,[r['pressure_max_bar'] for r in rows],label='Maximum');axes[2].set_ylabel('Pressure (bar absolute)');axes[2].legend()
    for ax in axes:ax.set_xlabel('Elapsed time (h)');ax.grid(alpha=.2)
    fig.suptitle('Assumed packed bed | 50°C, 120 bar, 4.5 molal NaCl | conversion half-time is a scenario')
    fig.savefig(output/'history.png',dpi=160);plt.close(fig)
    fig=plt.figure(figsize=(11,5),layout='constrained')
    for panel,(values,label) in enumerate([(concentration*1000,'H₂ concentration (g/m³ pore fluid)'),(m.reservoir_state()[:,0],'Pressure (bar absolute)')],1):
        ax=fig.add_subplot(1,2,panel,projection='3d');dots=ax.scatter(*xyz.T,c=values,cmap='viridis',s=85,depthshade=False)
        fig.colorbar(dots,ax=ax,shrink=.65,pad=.07,label=label)
        speed=np.linalg.norm(velocity,axis=1);scale=.25/max(speed.max(),1e-12)
        ax.quiver(*xyz.T,*(velocity*scale).T,color='#27333e',arrow_length_ratio=.25,linewidth=.8)
        for idx,name,color in [(0,'INJ','#1568a8'),(m.nb-1,'PROD','#b44232')]:
            ax.scatter(*xyz[idx],marker='*',s=190,color=color);ax.text(*xyz[idx],name)
        ax.set(xlim=(0,m.cfg['length_x_m']),ylim=(0,m.cfg['length_y_m']),zlim=(m.height,0),xlabel='x (m)',ylabel='y (m)',zlabel='Depth (m)')
        ax.set_xticks([0,1,2]);ax.set_yticks([0,1,2]);ax.set_zticks([0,.5,1.]);ax.tick_params(pad=1,labelsize=9)
        ax.xaxis.labelpad=5;ax.yaxis.labelpad=5
        ax.set_box_aspect((m.cfg['length_x_m'],m.cfg['length_y_m'],m.height));ax.view_init(elev=23,azim=-55)
    fig.suptitle(f'Native 3D packed-bed fields at {m.cfg["runtime_days"]*24:g} h; arrows show brine Darcy flow\nArrow lengths share one scale; maximum Darcy speed {speed.max():.3g} m/day')
    fig.savefig(output/'packed_bed_3d.png',dpi=180);plt.close(fig)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--suite',action='store_true');parser.add_argument('--output',type=Path);parser.add_argument('--config',type=Path)
    args=parser.parse_args();overrides=json.loads(args.config.read_text()) if args.config else {}
    root=args.output or Path(__file__).with_name('output')
    for name,extra in ([('reactive',{}),('zero_source',{'reaction_enabled':False})] if args.suite else [('reactive',{})]):
        summary,_=run_case(overrides|extra,root/name);print(json.dumps({k:summary[k] for k in ['generated_H2_kg','aqueous_H2_exported_kg','remaining_dissolved_H2_kg','max_native_relative_balance_error']},indent=2))

if __name__=='__main__':main()
