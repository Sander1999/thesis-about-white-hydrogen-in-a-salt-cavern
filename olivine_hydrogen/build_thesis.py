"""Synchronize calculated figures and compile the self-contained LaTeX source."""
import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent


def synchronize():
    source=HERE/'thesis';figures=source/'figures';figures.mkdir(exist_ok=True)
    for name in ['batch_scenarios.png','partition_recovery.png','energy_balance.png','sensitivity.png','alteration_feedbacks.png']:
        shutil.copy2(HERE/'results'/name,figures/name)
    for name in ['fe_reading_repeatability.pdf','fe_and_conditional_capacity.pdf']:
        shutil.copy2(HERE/'data/figures'/name,figures/name)
    shutil.copy2(HERE/'darts/output/reactive/packed_bed_3d.png',figures/'packed_bed_3d.png')
    for name in ['hydrogen_comparison.png','feedback_and_pumping.png']:
        shutil.copy2(HERE/'feedback/output'/name,figures/name)
    return source


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=HERE.parent/'bachelor_thesis_Sander_Bert_da_Costa.pdf')
    parser.add_argument('--compiler',default='/opt/anaconda3/bin/tectonic')
    parser.add_argument('--sync-only',action='store_true')
    args=parser.parse_args();source=synchronize()
    if args.sync_only:return
    with tempfile.TemporaryDirectory(prefix='olivine-thesis-build-') as temporary:
        subprocess.run([args.compiler,'--untrusted','--outdir',temporary,str(source/'main.tex')],check=True)
        args.output.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(Path(temporary)/'main.pdf',args.output)
    print(args.output)


if __name__=='__main__':main()
