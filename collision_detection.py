# Fernando Hueso-Gonzalez
# Massachusetts General Hospital and Harvard Medical School
#
# Illustrates how you can load a 3D model of the LINAC into RayStation for collision detection

from connect import *

fileName = '\\\\Client\\D$\\STL\\RotatingHeads.stl'

case = get_current('Case')
examinationName = 'CT 1'

# create ROI
roiName = 'ImportedModel2'
roiColor = 'Blue'
roiType = 'Support'
case.PatientModel.CreateRoi(Name=roiName, Color=roiColor, Type=roiType)

# import mesh from file
geo = case.PatientModel.StructureSets[examinationName].RoiGeometries[roiName]
geo.ImportRoiGeometryFromSTL(FileName=fileName, UnitInFile='Millimeter', TransformationMatrix=None)
