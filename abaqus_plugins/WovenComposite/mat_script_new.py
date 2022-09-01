import numpy as np
import Eshelby
from abaqus import *
from abaqusConstants import *

def assign_woven_composite_mat(modelName,partName,matrix_E,matrix_poissons,matrix_alpha,matrix_num_KB,matrix_cpt,matrix_tinit,matrix_viscoelastic,matrix_viscoplstic_damage,yarn_R22,yarn_R12,yarn_S,yarn_beta,yarn_gamma_inf,yarn_a22,yarn_a12,yarn_m_gauss,yarn_n_gauss,yarn_inclusion_dimension,yarn_stiffness):
	#
	#Matrix elements, yarn elements
	#
	matrix_elements = [ e.label for e in mdb.models[modelName].rootAssembly.sets['MATRIX'].elements]
	all_elements 	= [ e.label for e in mdb.models[modelName].rootAssembly.sets['ALL'].elements]
	yarn_elements 	= [ ele for ele in all_elements if ele not in matrix_elements]
	#
	#Model and part
	#
	model	= mdb.models[modelName]
	part	= model.parts[partName]
	#
	#Fetch discrete filed name which contains the material orientations exported from TexGen
	#
	fieldName = model.discreteFields.keys()[0]
	#
	#Fetch distribution data from discrete field
	#
	data 	= model.discreteFields[fieldName].data[0]
	domain 	= data.domain
	table 	= data.table
	width 	= data.dataWidth
	#
	#Yarn material constants
	#
	material_constants	= [yarn_R22,yarn_R12,yarn_S,yarn_beta,yarn_gamma_inf,yarn_a22,yarn_a12]
	initial_stiffness 	= [ele for row in yarn_stiffness for ele in row] # C1111,C1122,C2222,C2233,C1212
	C_0 				= np.zeros((6,6),dtype=float,order='F')
	#
	#row 1
	#
	C_0[0][0] = initial_stiffness[0]
	C_0[0][1] = initial_stiffness[1]
	C_0[0][2] = C_0[0][1]
	#
	#row 2
	#
	C_0[1][0] = C_0[0][1]
	C_0[1][1] = initial_stiffness[2]
	C_0[1][2] = initial_stiffness[3]
	#
	#row 3
	#
	C_0[2][0] = C_0[0][2]
	C_0[2][1] = C_0[1][2]
	C_0[2][2] = C_0[1][1]
	#
	#Lower diagonal
	#
	C_0[3][3] = initial_stiffness[4]
	C_0[4][4] = C_0[3][3]
	C_0[5][5] = 0.5 * (C_0[2][2] - C_0[1][2])
	#
	#Void inclusion dimension
	#
	a	= [ele for row in yarn_inclusion_dimension for ele in row]
	#
	#Compute eshelby tensor
	#
	SE	= Eshelby.eshelby_tensor(C_0,a,yarn_m_gauss,yarn_n_gauss)
	print('Eshelby tensor computed.')
	Eshelby_t = [SE[0][0],SE[0][1],SE[0][2],SE[1][0],SE[1][1],SE[1][2],SE[2][0],SE[2][1],SE[2][2],SE[3][3],SE[4][4],SE[5][5]]
	#
	#Loop over yarn elements and create material definition for each yarn element
	#
	for ele in yarn_elements:
		element_number	= ele
		set_name 		= 'Set_yarn_e_'+str(element_number)
		material_name 	= 'Mat_yarn_e_'+str(element_number)
		section_name 	= 'Section_yarn_e_'+str(element_number)		
		#
		#Create set
		#
		part.SetFromElementLabels(name=set_name,elementLabels=(element_number,))	
		#
		#Get material orientation
		#
		orientation = table[width*(element_number-1):width*element_number]
		#
		#Assemble material constants
		#
		constants =[2.0]
		constants.extend(orientation[:])
		constants.extend(initial_stiffness[:])
		constants.extend(Eshelby_t[:])
		constants.extend(material_constants[:])
		#
		#Create a material for yarn
		# 
		model.Material(name=material_name)
		model.materials[material_name].UserMaterial(mechanicalConstants=tuple(constants))
		model.materials[material_name].Depvar(n=14)
		#
		#Create section
		#
		model.HomogeneousSolidSection(name=section_name, material=material_name, thickness=None)
		#
		#Section assignment
		#
		region = part.sets[set_name]
		part.SectionAssignment(region=region, sectionName=section_name, offset=0.0,offsetType=MIDDLE_SURFACE, offsetField='',thicknessAssignment=FROM_SECTION)
	
	print('Yarn material assigned.')
	#
	#Define material and section names
	# 
	matrix_material_name	= 'mat_matrix'
	matrix_section_name		= 'Section_matrix'
	matrix_set_name			= 'Set_matrix'
	# 
	# Assemble material constants
	# 
	matrix_mat_constants		= [1.0,matrix_E,matrix_poissons,matrix_alpha,matrix_num_KB]
	matrix_thermal_const		= [0.0,0.0,0.0,0.0,0.0,0.0,matrix_cpt,matrix_tinit]
	kelvin_branch_param			= [ele for row in matrix_viscoelastic for ele in row]
	visco_plastic_damage_param	= [ele for row in matrix_viscoplstic_damage for ele in row]
	matrix_mat_constants.extend(kelvin_branch_param[:])
	matrix_mat_constants.extend(visco_plastic_damage_param[:])
	matrix_mat_constants.extend(matrix_thermal_const[:])
	# 
	#Create material
	# 
	model.Material(name=matrix_material_name)
	model.materials[matrix_material_name].UserMaterial(mechanicalConstants=tuple(matrix_mat_constants))
	model.materials[matrix_material_name].Depvar(n=36)
	#
	#Create set for matrix
	#
	model	= mdb.models[modelName]
	part	= model.parts[partName]
	part.SetFromElementLabels(name=matrix_set_name,elementLabels=tuple(matrix_elements))
	#
	#Create section
	#
	model.HomogeneousSolidSection(name=matrix_section_name, material=matrix_material_name, thickness=None)
	#
	#Section assignment
	#
	region = mdb.models[modelName].parts[partName].sets[matrix_set_name]
	part.SectionAssignment(region=region, sectionName=matrix_section_name, offset=0.0,offsetType=MIDDLE_SURFACE, offsetField='',thicknessAssignment=FROM_SECTION)
	print('Matrix material assigned.')
	print('Woven composite material assignment completed.')

