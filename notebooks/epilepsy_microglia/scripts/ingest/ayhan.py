"""Offline Ayhan ingestion from GEO CSV and SOFT metadata."""
import argparse,gzip,json,sys,csv
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from epilepsy_microglia.io import write_study
from epilepsy_microglia.metadata import standardize_obs
from epilepsy_microglia.validation import validate_adata
from epilepsy_microglia.download import sha256
def parse_geo(path):
 with gzip.open(path,'rt') as f:t=f.read()
 rows=[]
 for b in t.split('^SAMPLE = ')[1:]:
  r={'source_sample_id':b.splitlines()[0]}
  for l in b.splitlines():
   if l.startswith('!Sample_title = '):r['sample_id']=l.split(' = ',1)[1]
   elif l.startswith('!Sample_source_name_ch1 = '):r['donor_id']=l.split(' = ',1)[1]
   elif l.startswith('!Sample_characteristics_ch1 = '):
    k,v=l.split(' = ',1)[1].split(': ',1);r[k.replace(' ','_').lower()]=v
  rows.append(r)
 return pd.DataFrame(rows)
def run(root):
 raw=root/'data/raw/ayhan';meta=root/'data/metadata/ayhan';out=root/'data/processed/ayhan';inv=json.loads((meta/'download_inventory.json').read_text())
 for z in inv:
  if sha256(raw/z['filename'])!=z['sha256']:raise ValueError('raw checksum changed')
 g=parse_geo(raw/'GSE160189_family.soft.gz');cs=[z for z in inv if z['filename'].endswith(('.csv','.csv.gz'))]
 if len(cs)!=1:raise ValueError('expected one count CSV')
 with gzip.open(raw/cs[0]['filename'],'rt',newline='') as h:
  reader=csv.reader(h); header=next(reader); cells=pd.Index(header[1:],name='source_cell_id'); rows=[]; cols=[]; vals=[]; names=[]
  for i,row in enumerate(reader):
   names.append(row[0])
   for j,v in enumerate(row[1:]):
    if v and v!='0': rows.append(j); cols.append(i); vals.append(int(v))
 x=sparse.coo_matrix((vals,(rows,cols)),shape=(len(cells),len(names)),dtype=np.int32).tocsr();genes=pd.Series(names,dtype=str)
 var=pd.DataFrame(index=pd.Index(genes,name='gene_name'));var['source_gene_id']=genes.to_numpy();obs=standardize_obs(pd.DataFrame(index=cells),'ayhan')
 for c in g:obs[c]=pd.NA
 obs['source_file']='ayhan/'+cs[0]['filename'];obs['source_accession']='GSE160189';obs['sample_mapping_available']=False
 a=ad.AnnData(X=x,obs=obs,var=var);a.layers['counts']=a.X;a.uns['ingestion']={'raw_counts_available':True,'expression_type':'raw_counts','source_files':[z['filename'] for z in inv],'metadata_sources':['GSE160189_family.soft.gz'],'transformations':'CSV transpose only; no QC or normalization'}
 out.mkdir(parents=True,exist_ok=True);g.to_csv(meta/'sample_metadata.csv',index=False);dest=write_study(a,study='ayhan',project_root=root);r=ad.read_h5ad(dest);e=validate_adata(r,study='ayhan')
 if e:raise ValueError('\n'.join(e))
 print(dest,r.shape,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,default=ROOT);run(p.parse_args().project_root.resolve())
