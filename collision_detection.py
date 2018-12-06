# Erik's example script
#
# Illustrates how you can load STL geometries into RayStation

from connect import *

fileName = 'C:\\file.STL'

case = get_current('Case')
examinationName = 'CT 1'

# create ROI
roiName = 'ImportedModel'
roiColor = 'Red'
roiType = 'Support'
case.PatientModel.CreateRoi(Name=roiName, Color=roiColor, Type=roiType)

# import mesh from file
geo = case.PatientModel.StructureSets[examinationName].RoiGeometries[roiName]
geo.ImportRoiGeometryFromSTL(FileName=fileName, UnitInFile='Millimeter', TransformationMatrix=None)
