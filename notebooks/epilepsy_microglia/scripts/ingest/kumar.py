"""Offline Kumar ingestion from original GEO MTX/TSV files and SOFT metadata."""
import argparse,gzip,json,sys
from pathlib import Path
import anndata as ad
import pandas as pd
from scipy import sparse
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from epilepsy_microglia.io import write_study
from epilepsy_microglia.metadata import standardize_obs
from epilepsy_microglia.validation import validate_adata
from epilepsy_microglia.kumar import read_triplet
def parse_geo(path):
 with gzip.open(path,'rt') as f:t=f.read()
 rows=[]
 for b in t.split('^SAMPLE = ')[1:]:
  r={'source_sample_id':b.splitlines()[0]}
  for l in b.splitlines():
   if l.startswith('!Sample_title = '):r['sample_id']=l.split(' = ',1)[1]
   elif l.startswith('!Sample_source_name_ch1 = '):r['tissue']=l.split(' = ',1)[1]
  rows.append(r)
 return pd.DataFrame(rows)
def run(root):
 raw=root/'data/raw/kumar';meta=root/'data/metadata/kumar';out=root/'data/processed/kumar';inv=json.loads((meta/'download_inventory.json').read_text());g=parse_geo(raw/'GSE201048_family.soft.gz')
 if len(g)!=11:raise ValueError('GEO does not contain 11 samples')
 files={r['source_sample_id']:[] for r in g.to_dict('records')}
 for z in inv:
  if z['filename'].endswith(('.mtx.gz','.tsv.gz')):
   gsm=z['filename'].split('_',1)[0]; files.setdefault(gsm,[]).append(z['filename'])
 mats=[];obsall=[];var0=None
 for r in g.to_dict('records'):
  full,bar,feat=read_triplet(raw,files[r['source_sample_id']]);mask=feat.feature_type.eq('Gene Expression').to_numpy();genes=feat.loc[mask].copy();genes.index=pd.Index(genes.source_gene_id,name='gene_id')
  if var0 is None:var0=genes
  elif not genes.equals(var0):raise ValueError('gene axes differ')
  o=standardize_obs(pd.DataFrame(index=bar),'kumar');o['sample_id']=r['sample_id'];o['source_sample_id']=r['source_sample_id'];o['tissue']=r['tissue'];o['donor_id']=pd.NA;o['diagnosis']=pd.NA;o['source_file']='kumar/'+r['source_sample_id'];o['assay']='scRNA-seq';o['platform']='10x Genomics Chromium + CITE-seq';o['author_annotation_available']=False;o.index=pd.Index([r['source_sample_id']+':'+x for x in bar],name='cell_id');mats.append(full[:,mask].tocsr());obsall.append(o)
 a=ad.AnnData(X=sparse.vstack(mats,format='csr'),obs=pd.concat(obsall),var=var0);a.layers['counts']=a.X;a.uns['ingestion']={'raw_counts_available':True,'expression_type':'raw_counts','source_files':[z['filename'] for z in inv],'metadata_sources':['GSE201048_family.soft.gz'],'transformations':'MTX transpose and RNA feature selection only; no QC or normalization'}
 out.mkdir(parents=True,exist_ok=True);dest=write_study(a,study='kumar',project_root=root);r=ad.read_h5ad(dest);e=validate_adata(r,study='kumar')
 if e:raise ValueError('\n'.join(e))
 print(dest,r.shape,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,default=ROOT);run(p.parse_args().project_root.resolve())
