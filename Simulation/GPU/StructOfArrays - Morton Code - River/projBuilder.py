

import math

from pathlib import Path
from string import Template


def getFileString(path):
	return	Path(path).read_text(encoding="utf8")

def buildVoxelProj():
	print('Reading Voxel Proj')
	doc = getFileString('templates/voxelproj.xml')

	numParticles = 64000
	voxeldims    = (290,268,386)
	#voxeldims    = (0.7*290,0.7*268,0.7*386)
	voxeldims    = (290*1.5,268*1.5,386*1.5)
	
	variables =  {
		'NUMPARTICLES': numParticles   ,
		'NUMVERTS'    : numParticles*6 ,
		'WORKGROUPS'  : math.ceil(numParticles/32),
		'VOX_SIDE_X'  : voxeldims[0],
		'VOX_SIDE_Y'  : voxeldims[1],
		'VOX_SIDE_Z'  : voxeldims[2],
		'VOX_SIDE_GROUPS_X' : math.ceil(voxeldims[0]/32),
		'VOX_SIDE_GROUPS_Y' : math.ceil(voxeldims[1]/32),
		'windowSimul'   :   getFileString('templates\simulationInterface.xml'),
		'simpipeline'	:	getFileString('templates\simpipeline.xml'),
		'simattributes'	:	getFileString('templates\simattributes.xml'),
	}

	t = Template(doc)
	outpath = 'voxel.xml'
	print(f'Writing Voxel Proj to {outpath}')
	Path(outpath).write_text(t.substitute(variables),encoding='utf8')


def buildOldProj():
	print('Reading Old Proj templates/oldRender.xml')
	doc = getFileString('templates/oldRender.xml')


	numParticles = 64000
	voxeldims    = (290,268,386)
	
	variables =  {
		'NUMPARTICLES': numParticles   ,
		'NUMVERTS'    : numParticles*6 ,
		'WORKGROUPS'  : math.ceil(numParticles/32),
		'VOX_SIDE_X'  : voxeldims[0],
		'VOX_SIDE_Y'  : voxeldims[1],
		'VOX_SIDE_Z'  : voxeldims[2],
		'VOX_SIDE_GROUPS_X' : math.ceil(voxeldims[0]/32),
		'VOX_SIDE_GROUPS_Y' : math.ceil(voxeldims[1]/32),
		'windowSimul'   :   getFileString('templates\simulationInterface.xml'),
		'simpipeline'	:	getFileString('templates\simpipeline.xml'),
		'simattributes'	:	getFileString('templates\simattributes.xml'),
	}

	t = Template(doc)
	outpath = 'proj.xml'
	print(f'Writing Old Proj to {outpath}')
	Path(outpath).write_text(t.substitute(variables),encoding='utf8')

buildOldProj()
buildVoxelProj()