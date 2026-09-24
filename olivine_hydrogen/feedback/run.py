"""Run stationary/recirculating/renewing reduced-feedback comparisons."""
from dataclasses import replace
from pathlib import Path
import argparse
import csv
import json
import numpy as np
from .model import Config, MODES, simulate


def scenarios(cfg):
    return {
        "no_feedback_null": replace(cfg,solid_expansion_fraction=0.,maximum_area_factor=1.,
                                    passivation_coefficient=0.,porosity_access_exponent=0.),
        "intrinsic_limited": cfg,
        "film_only_hypothesis": replace(cfg,film_damkohler=1.),
        "inhibition_only_hypothesis": replace(cfg,product_inhibition_coefficient=1.),
        "contact_inhibition_hypothesis": replace(cfg,film_damkohler=1.,product_inhibition_coefficient=1.),
        "closure_hypothesis": replace(cfg,film_damkohler=1.,product_inhibition_coefficient=1.,
                                      maximum_area_factor=1.,crack_relief_fraction=0.,passivation_coefficient=4.),
        "opening_hypothesis": replace(cfg,film_damkohler=1.,product_inhibition_coefficient=1.,
                                      crack_relief_fraction=.75,passivation_coefficient=.5),
        "zero_source": replace(cfg,rate_multiplier=0.),
    }


def _write_csv(path, rows):
    with path.open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def make_plots(histories, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors={"stationary":"#397da2","recirculating":"#62856c","once_through":"#b36e38"}
    with plt.rc_context({"font.size":16,"axes.labelsize":14,"axes.titlesize":15,
                        "xtick.labelsize":14,"ytick.labelsize":14,
                        "axes.spines.top":False,"axes.spines.right":False}):
        fig,axes=plt.subplots(2,3,figsize=(12,8.5),layout="constrained")
        for row,scenario in enumerate(("intrinsic_limited","contact_inhibition_hypothesis")):
            for mode in MODES:
                h=histories[(scenario,mode)]
                for col,(key,title) in enumerate((("generated_h2_kg","Generated hydrogen"),("retained_h2_kg","Retained in reactor"),("total_deliverable_h2_kg","Collectable hydrogen\nOutlet + terminal flash"))):
                    axes[row,col].plot(h["time_days"],h[key],label=mode.replace("_"," "),color=colors[mode],lw=2,
                                       ls={"stationary":"-","recirculating":"--","once_through":":"}[mode])
                    axes[row,col].set(xlabel="Time (days)",ylabel="Hydrogen (kg)",title=title);axes[row,col].grid(alpha=.2)
            axes[row,0].text(.02,.98,"Intrinsic rate limit" if row==0 else "Contact + inhibition",
                             transform=axes[row,0].transAxes,va="top",fontsize=14)
        axes[0,2].legend(loc="upper left",fontsize=14)
        fig.suptitle("Flow, reaction and hydrogen recovery",fontsize=18)
        fig.savefig(output/"hydrogen_comparison.png",dpi=180);plt.close(fig)
        fig,axes=plt.subplots(1,3,figsize=(12,4.8),layout="constrained")
        for mode in MODES:
            h=histories[("contact_inhibition_hypothesis",mode)]
            for ax,key,title in zip(axes,("porosity","area_factor","flow_pump_kwh"),("Reactive-bed\nporosity","Accessible surface\n(relative area)","Flow pumping\n(cumulative kWh)")):
                ax.plot(h["time_days"],h[key],label=mode.replace("_"," "),color=colors[mode],lw=2)
                ax.set(xlabel="Time (days)",title=title);ax.grid(alpha=.2)
        axes[0].legend(fontsize=14);fig.suptitle("Assumed expansion and cracking feedbacks",fontsize=18)
        fig.savefig(output/"feedback_and_pumping.png",dpi=180);plt.close(fig)


def run_suite(output=None,cfg=Config(),plots=True):
    output=Path(output or Path(__file__).with_name("output"));output.mkdir(parents=True,exist_ok=True)
    histories={};summaries={};endpoints=[]
    for scenario,scenario_cfg in scenarios(cfg).items():
        for mode in MODES:
            h,s=simulate(mode,scenario_cfg);histories[(scenario,mode)]=h
            summaries[f"{scenario}/{mode}"]=s
            rows=[{"scenario":scenario,"mode":mode,**{k:float(v[i]) for k,v in h.items()}} for i in range(len(h["time_days"]))]
            _write_csv(output/f"{scenario}_{mode}.csv",rows)
            for t in (3.,30.,90.,365.,cfg.duration_days):
                idx=np.flatnonzero(np.isclose(h["time_days"],t))
                if len(idx) and not any(r["scenario"]==scenario and r["mode"]==mode and r["time_days"]==t for r in endpoints):
                    endpoints.append(rows[int(idx[0])])
    _write_csv(output/"comparison_endpoints.csv",endpoints)
    ratios={}
    for scenario in scenarios(cfg):
        ref=summaries[f"{scenario}/stationary"]["final"]
        ratios[scenario]={}
        for mode in MODES:
            last=summaries[f"{scenario}/{mode}"]["final"]
            ratios[scenario][mode]={"generated_vs_stationary":last["generated_h2_kg"]/ref["generated_h2_kg"] if ref["generated_h2_kg"] else None,
                                   "retained_vs_stationary":last["retained_h2_kg"]/ref["retained_h2_kg"] if ref["retained_h2_kg"] else None}
    result={"cases":summaries,"end_time_ratios":ratios,
            "warning":"All feedbacks/rates are uncalibrated scenarios; apparent rankings depend on the assumptions. Terminal flashes are independent endpoints, not cumulative repeated recovery."}
    (output/"summary.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    if plots:make_plots(histories,output)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--days",type=float,default=365.)
    parser.add_argument("--no-plots",action="store_true")
    args=parser.parse_args()
    result=run_suite(args.output,replace(Config(),duration_days=args.days),not args.no_plots)
    for key,value in result["cases"].items():
        f=value["final"]
        print(f"{key}: generated={f['generated_h2_kg']:.6g} kg, retained={f['retained_h2_kg']:.6g} kg, collectable={f['total_deliverable_h2_kg']:.6g} kg; additional water={f['feed_water_kg']:.6g} kg")
