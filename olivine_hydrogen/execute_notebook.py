"""Execute the notebook with this exact interpreter using a temporary kernel."""
import json
import os
import sys
import tempfile
from pathlib import Path
import nbformat
from nbclient import NotebookClient


def main():
    here=Path(__file__).resolve().parent
    path=here/'Olivine_data_transformation.ipynb'
    notebook=nbformat.read(path,as_version=4)
    previous=os.environ.get('JUPYTER_PATH')
    with tempfile.TemporaryDirectory(prefix='olivine-kernel-') as directory:
        spec=Path(directory)/'kernels'/'olivine-darts';spec.mkdir(parents=True)
        (spec/'kernel.json').write_text(json.dumps({
            'argv':[sys.executable,'-B','-m','ipykernel_launcher','-f','{connection_file}'],
            'display_name':'Olivine model (shared DARTS Python)', 'language':'python'}))
        os.environ['JUPYTER_PATH']=directory+(os.pathsep+previous if previous else '')
        try:
            NotebookClient(notebook,kernel_name='olivine-darts',timeout=600,
                resources={'metadata':{'path':str(here.parent)}}).execute()
            errors=[o for cell in notebook.cells if cell.cell_type=='code'
                    for o in cell.get('outputs',[]) if o.output_type=='error']
            if errors: raise RuntimeError('Notebook contains execution errors')
            nbformat.write(notebook,path)
        finally:
            if previous is None:os.environ.pop('JUPYTER_PATH',None)
            else:os.environ['JUPYTER_PATH']=previous
    print('Executed:',path)


if __name__=='__main__':main()
