
import pandas as pd
import requests
import glob,os

swc_url_base='https://celltypes.brain-map.org/api/v2/well_known_file_download/{}'


def allen_retrieve_swc(fCells_path,SWC_path):
	#use a table of a structure of the
	#table provided by Allen Brain
	#which includes mouse cell metadata
	#to retrieve dendrites structure
	#for all cells
	#input: filtered table path
	#       path for swc files
	#stores swc files with the names of
	#  %speciment_id%.swc in the swc path
	dt = pd.read_csv(fCells_path) 
	SWC_LUT=dt[['specimen__id','nrwkf__id']]
	
	swc_path_base=SWC_PATH+'/{}.swc'
	for i,lut in SWC_LUT.iterrows():#
		sp_id=lut['specimen__id']
		swc_fid=lut['nrwkf__id']
		swc_url=swc_url_base.format(swc_fid)
		swc_path=swc_path_base.format(sp_id)
		r = requests.get(swc_url)

		with open(swc_path, 'wb') as f:
			f.write(r.content)
   

#list_local_cells_allen()
