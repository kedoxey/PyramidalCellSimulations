import requests
import json
import os
from zipfile import ZipFile
from urllib.request import urlopen

cwd = os.getcwd()


def compile_mechs(cwd, hocs_dir, mod_dir):
    if not os.path.exists(os.path.join(hocs_dir, 'x86_64')):
        os.chdir(hocs_dir)

        os.system(f'nrnivmodl {mod_dir}')

        os.chdir(cwd)
    else:
        print('Mechanisms already compiled!')


# download specified version of model from neuroml-db 
def download_from_nmldb(model_id, version):
    
    zip_url = f'https://neuroml-db.org/GetModelZip?modelID={model_id}&version={version}'

    unzip_path = os.path.join(cwd,'models',version,f'{model_id}-{version}')
    zip_path = os.path.join(cwd,'models','zips',model_id+'.zip')

    if not os.path.exists(unzip_path):
        # download model if it has not been already
        nmldb_response = requests.get(zip_url)
        open(zip_path,'wb').write(nmldb_response.content)

        os.makedirs(unzip_path)

        with ZipFile(zip_path,'r') as zObject:
            zObject.extractall(path=unzip_path)

        print(f'Model {model_id} successfully downloaded!')
    else:
        print(f'Model {model_id} already downloaded.')

    
