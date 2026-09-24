"""Regenerate the batch model, numerical tables and publication figures."""
import csv
import json
from dataclasses import replace
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .batch import Setup, SCENARIOS, evaluate, metadata, capacity_mol, M_H2

HERE = Path(__file__).resolve().parent


def fresh_eifel_setup():
    with (HERE/'data/processed/subgroup_summary.csv').open() as f:
        rows=list(csv.DictReader(f))
    # Export subgroup labels are reconstructed from the original notebook.
    matches=[r for r in rows if r['sample']=='sample5' and r['subgroup']=='fresh surface']
    if len(matches)!=1:
        raise ValueError('Expected exactly one sample 5 fresh subgroup')
    return replace(Setup(), fe_mass_fraction=float(matches[0]['Fe_mean_ppm'])/1e6)


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def run(output=HERE/'results'):
    output=Path(output); output.mkdir(parents=True,exist_ok=True)
    cfg=fresh_eifel_setup()
    times=np.unique(np.r_[0.,3.,30.,365.,3650., np.geomspace(.05,3650,180)])
    histories=[]; endpoints=[]
    for name,half in SCENARIOS.items():
        history=evaluate(times,half,cfg)
        for i,t in enumerate(times):
            row={'scenario':name, **{k:float(v[i]) for k,v in history.items()}}
            histories.append(row)
            if t in (3,365,3650): endpoints.append(row)
    write_csv(output/'batch_history.csv',histories)
    write_csv(output/'batch_endpoints.csv',endpoints)
    ceiling={k:float(v[0]) for k,v in evaluate([1e9],3,cfg).items()}
    meta=metadata(cfg); meta['complete_conversion_endpoint']=ceiling
    meta['engineering_exclusions']=['heat losses and mineral reaction enthalpy',
        'rock mining, grinding, transport and replacement', 'gas drying and purification',
        'cavern construction, geomechanics, microbial losses and well lifting work',
        'additional fluid circulation and compression beyond specified boundary']
    (output/'batch_summary.json').write_text(json.dumps(meta,indent=2)+'\n')

    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'savefig.dpi':180,'figure.constrained_layout.use':True})
    colors=['#555555','#397ca8','#df8c35','#903c60']
    labels=['Zero production','Half-time 10 years','Half-time 1 year','Half-time 3 days (target)']
    fig,axs=plt.subplots(1,2,figsize=(11,4.3))
    for (name,half),color,label in zip(SCENARIOS.items(),colors,labels):
        x=evaluate(times,half,cfg)
        axs[0].plot(times,x['generated_h2_kg'],label=label,color=color)
        axs[1].plot(times,x['gross_electricity_kwh'],color=color,label=label)
    for ax in axs:
        ax.set_xscale('symlog',linthresh=3);ax.set_xlabel('Batch duration (days)');ax.grid(alpha=.2)
    axs[0].axhline(capacity_mol(cfg)*M_H2,ls=':',color='black',lw=1)
    axs[0].set_ylabel('Generated hydrogen (kg per 10 t rock)')
    axs[1].set_ylabel('Gross electricity after recovery (kWh)')
    axs[0].legend(fontsize=8,loc='lower right')
    fig.suptitle('Finite iron inventory and assumed reaction times')
    fig.savefig(output/'batch_scenarios.png');plt.close(fig)

    fig,axs=plt.subplots(1,2,figsize=(11,4.2))
    x=evaluate(times,365,cfg)
    axs[0].stackplot(times,x['dissolved_reactor_h2_kg'],x['free_reactor_h2_kg'],
                     labels=['Dissolved at 120 bar','Free gas at 120 bar'],colors=['#397ca8','#e4ae55'])
    axs[1].stackplot(times,x['retained_separator_h2_kg'],x['uncollected_separator_h2_kg'],x['collected_h2_kg'],
                     labels=['Dissolved at 10 bar','Uncollected gas','Collected gas'],colors=['#397ca8','#aaaaaa','#609276'])
    for ax in axs:
        ax.set_xlabel('Batch duration (days)');ax.set_ylabel('Hydrogen (kg)');ax.legend(fontsize=8,loc='upper left')
        ax.set_xlim(0,1500);ax.grid(alpha=.2)
    fig.suptitle('Where the hydrogen goes: assumed one-year reaction half-time')
    fig.savefig(output/'partition_recovery.png');plt.close(fig)

    names=['3-day target\nat 72 hours','1-year half-time\nat 1 year','Complete\nconversion ceiling']
    states=[evaluate([3],3,cfg),evaluate([365],365,cfg),evaluate([1e9],3,cfg)]
    fig,ax=plt.subplots(figsize=(9,4.5)); ids=np.arange(3);width=.24
    for delta,key,label,color in [(-width,'gross_electricity_kwh','Gross electricity','#609276'),
             (0,'warm_pressurized_balance_kwh','Warm, pressurized boundary balance','#397ca8'),
             (width,'standalone_balance_kwh','Electric heat + surface pressurization','#b55453')]:
        ax.bar(ids+delta,[r[key][0] for r in states],width,label=label,color=color)
    ax.axhline(0,color='black',lw=.7);ax.set_xticks(ids,names);ax.set_ylabel('Electricity per 10 t rock (kWh)')
    ax.legend(fontsize=8,loc='lower left');ax.grid(axis='y',alpha=.15)
    ax.set_title('Energy result depends on the operating boundary')
    fig.savefig(output/'energy_balance.png');plt.close(fig)

    # Accessible iron is a major uncertainty, rather than just XRF repeatability.
    fig,ax=plt.subplots(figsize=(8,4.3))
    fractions=np.linspace(0,1,101)
    for m,color in zip([3.5,4.5,4.95],['#397ca8','#609276','#df8c35']):
        y=[evaluate([365],365,replace(cfg,accessible_fraction=f,salt_molality=m))['gross_electricity_kwh'][0] for f in fractions]
        ax.plot(fractions,y,label=f'{m:g} mol NaCl/kg water',color=color)
    ax.set(xlabel='Fraction of iron accessible (Fe(II) fraction assumed 1)',ylabel='Gross electricity after 1 year (kWh)',
           title='Sensitivity to accessible iron and brine salinity')
    ax.legend();ax.grid(alpha=.2);fig.savefig(output/'sensitivity.png');plt.close(fig)
    from .plot_alteration import plot as plot_alteration
    plot_alteration(output)
    print(json.dumps({'ceiling_h2_kg':meta['magnetite_h2_ceiling_kg'],
                      'ceiling_gross_kwh':ceiling['gross_electricity_kwh'],
                      'heat_kwh_thermal':ceiling['sensible_heat_kwh_thermal'],
                      'pump_kwh':ceiling['initial_pump_electricity_kwh']},indent=2))
    return meta,endpoints


if __name__=='__main__': run()
