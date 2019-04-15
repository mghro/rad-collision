# collision_detection.py
#
# Fernando Hueso-Gonzalez - fhuesogonzalez@mgh.harvard.edu
# Massachusetts General Hospital and Harvard Medical School
#
# Loads a 3D model of a LINAC into RayStation for collision detection
#
# DICOM (RayStation) coordinates are as follows:
# x: from Right to Left
# y: from Anterior to Posterior
# z: from Inferior to Superior
#
# This is independent on the patient orientation (HFS, FFS, HFP, FFP).
# See http://dosimetry.dotdecimal.com/doku.php?id=dosimetry:userguide:proton_delivery_system_conventions
#
# In the transverse 2D viewer of RayStation, for FFS patient orientation, the coordinate system looks as follows:
# x <------|
#          |
#          |
#          V
#          y
# and z completes the triad (pointing to you, reader)
# However, it looks different for HFS orientation or any other, as they maintain the CT acquisition orientation
# and move instead the labels R, L, A, P, S, I.
#
# Let us define some rotation angles:
# a: rotation around z axis (from x to y axis), represented by rotation matrix R_z, with R_z[row=1,column=2] = -sin(a)
# b: rotation around y axis (from x to z axis), represented by rotation matrix R_y, with R_z[row=1,column=3] = -sin(b)
#
# g: DICOM gantry rotation (when couch angle is zero) is from y to x axis (opposite sign than a)
# c: Couch support rotation (when gantry angle is zero) is from x to z axis (same sign than b)
#
# 3D model coordinates
# It is assumed that the origin of the STL 3D models is exactly at room isocenter.
# Also, the orientation should be such that, when opening the STL file with Meshlab, you are looking towards the LINAC
# as if you were standing in the treatment room in front of it (no couch support rotation),
# and the gantry is at zero degrees.
# You might need to cleanup your model with Meshlab via opening it, click on
# Unify duplicated Vertices, Ok, File, Export Mesh As, .stl
# in order to avoid an error when importing in RayStation
#
# The initial affine transformation applied will be to rotate LINAC model according to gantry angle,
# then simulate couch angle via a negative couch rotation of the model,
# and finally translation of the model to the isocenter (iso.x,iso.y,iso.z) in the CT patient
# resulting in TransformationMatrix M(iso.x,iso.y,iso.z,c,g) = T(iso.x,iso.y,iso.z) * R_y ( b = -c ) * R_z ( a = -g ) =
# {'M11':cos(a)*cos(b), 'M12':-sin(a)*cos(b), 'M13':-sin(b), 'M14':iso.x,
#  'M21':sin(a)       , 'M22': cos(a)       , 'M23': 0     , 'M24':iso.y,
#  'M31':cos(a)*sin(b), 'M32':-sin(a)*sin(b), 'M33': cos(b), 'M34':iso.z,
#  'M41':0            , 'M42':0             , 'M43': 0     , 'M44':1          }
#
# See https://math.stackexchange.com/questions/2093314/rotation-matrix-of-rotation-around-a-point-other-than-the-origin
# See http://www2.clarku.edu/faculty/djoyce/trig/identities.html
#
# Subsequent input angles (g2,c2) will be computed as a differential affine transformation on top of the previous one (g,c),
# by translating the model isocenter back to the coordinate system origin,
# undoing the previous couch rotation (c)
# undoing the previous gantry rotation (g)
# and computing the previous TransformationMatrix M with the new angles. (g2,c2).
# Altogether, it yields the TransformationMatrix D = M(iso.x,iso.y,iso.z,c2,g2) * inverse(M(iso.x,iso.y,iso.z,c2,g2)) =
# M(iso.x,iso.y,iso.z,c2,g2) * R_z ( a = g ) * R_y ( b = c ) * T(-iso.x,-iso.y,-iso.z) =
# T(iso.x,iso.y,iso.z) * R_y ( b2 = -c2 ) * R_z ( a2 = -g2)  * R_z ( a = g ) * R_y ( b = c ) * T(-iso.x,-iso.y,-iso.z) =
# T(iso.x,iso.y,iso.z) * R_y ( b2 = -c2 ) * R_z ( d = g-g2) * R_y ( b = c ) * T(-iso.x,-iso.y,-iso.z) =
# resulting in TransformationMatrix D(iso.x,iso.y,iso.z,c2,g2,c,g) =
# {'M11':cos(d)*cos(b)*cos(b2)-sin(b)*sin(b2), 'M12':-sin(d)*cos(b2), 'M13':-cos(d)*sin(b)*cos(b2)-cos(b)*sin(b2), 'M14':iso.x-iso.x*(cos(d)*cos(b)*cos(b2)-sin(b)*sin(b2))+iso.y*sin(d)*cos(b2)+iso.z*(cos(d)*sin(b)*cos(b2)+cos(b)*sin(b2)),
#  'M21':sin(d)*cos(b)                       , 'M22': cos(d)        , 'M23':-sin(d)*sin(b)                       , 'M24':iso.y-iso.x* sin(d)*cos(b)                        -iso.y*cos(d)        +iso.z* sin(d)*sin(b)                        ,
#  'M31':cos(d)*cos(b)*sin(b2)+sin(b)*cos(b2), 'M32':-sin(d)*sin(b2), 'M33':-cos(d)*sin(b)*sin(b2)+cos(b)*cos(b2), 'M34':iso.z-iso.x*(cos(d)*cos(b)*sin(b2)+sin(b)*cos(b2))+iso.y*sin(d)*sin(b2)+iso.z*(cos(d)*sin(b)*sin(b2)-cos(b)*cos(b2)),
#  'M41':0                                   , 'M42': 0             , 'M43': 0                                   , 'M44':1                                                             }

class Part:

  def __init__(self, name, file, color):
    self.name = name
    self.file = file
    self.color=color

class Machine:

  def __init__(self, name, path, parts):
    self.name = name
    self.path = path
    self.parts = parts

#agility = Machine("Agility", "\\\\Client\\D$\\Dropbox (Partners Healthcare)\\LINAC\\Elekta Agility\\STL parts\\",[Part("Nozzle","RotatingHeads.stl","Blue") ])
#agility = Machine("Agility", "\\\\Client\\C$\\Users\\fh969\\Desktop\\STL parts\\",[Part("Nozzle","RotatingHeads.stl","Blue") ])
#agility = Machine("Agility", "\\\\Client\\H$\\STL parts\\",[Part("Nozzle","RotatingHeads.stl","Blue") ])
agility = Machine("Elekta Agility"   , "Y:\\STL parts\\Elekta Agility\\"   ,[Part("Nozzle","RotatingHeads.stl","Blue") ])
trilogy = Machine("Varian IX Trilogy", "Y:\\STL parts\\Varian IX Trilogy\\",[Part("Nozzle","RotatingHeads.stl","Blue") ])
truebeam= Machine("Varian TrueBeam"  , "Y:\\STL parts\\Varian TrueBeam\\"  ,[Part("Nozzle","RotatingHeads.stl","Blue") ])

machines = {agility.name : agility, trilogy.name : trilogy, truebeam.name: truebeam}

from math import cos, sin, radians, degrees, hypot, atan2
import sys
import ScriptClient

from connect import *

import clr
clr.AddReference("System.Windows.Forms")
clr.AddReference("System.Drawing")

from System.Windows.Forms import Application, Form, Label, ComboBox, Button, TextBox, TrackBar, FormStartPosition, TickStyle, Keys
from System.Drawing import Point, Size


# Define Forms class that will prompt the user to select a
# machine for creating a 3D model around isocenter
class SelectMachineForm(Form):

  def __init__(self):
    # Set the size of the form
    self.Size = Size(300, 200)
    # Set title of the form
    self.Text = 'Select Machine'

    # Add a label
    label = Label()
    label.Text = 'Please select a LINAC model'
    label.Location = Point(15, 15)
    label.AutoSize = True
    self.Controls.Add(label)

    # Add a ComboBox that will display the Machines to select
    # Define the items to show up in the combobox
    self.combobox = ComboBox()
    self.combobox.DataSource = machines.keys()
    self.combobox.Location = Point(15, 60)
    self.combobox.AutoSize = True
    self.Controls.Add(self.combobox)

    # Add button to press OK and close the form
    button = Button()
    button.Text = 'OK'
    button.AutoSize = True
    button.Location = Point(15, 100)
    button.Click += self.ok_button_clicked
    self.Controls.Add(button)

  def ok_button_clicked(self, sender, event):
    # Method invoked when the button is clicked
    # Save the selected machine name
    self.mach_name = self.combobox.SelectedValue
    # Close the form
    self.Close()

class SelectAngleForm(Form):
  def __init__(self):
    # Set the size of the form
    self.StartupPosition = FormStartPosition.Manual
    self.Location = Point(500,15)
    self.Size = Size(500, 300)
    # Set title of the form
    self.Text = 'Select Beam and Couch Angles'
    self.TopMost = True

    # Add a beam label
    labelB = Label()
    labelB.Text = 'Please select a beam angle in DEG [0:360].'
    labelB.Location = Point(15, 15)
    labelB.AutoSize = True
    self.Controls.Add(labelB)

    # Add a text box that to write the desired beam angle
    self.tboxB = TextBox()
    self.tboxB.Location = Point(15, 60)
    self.tboxB.Width = 55
    self.tboxB.Text = "0"
    self.tboxB.KeyDown += self.on_enter
    self.Controls.Add(self.tboxB)

    # Add a trackbar to slide to the desired beam angle
    self.tbB = TrackBar()
    self.tbB.TickStyle = TickStyle.Both
    self.tbB.TickFrequency = 10;
    self.tbB.Minimum = 0
    self.tbB.Maximum = 360
    self.tbB.Value = 0
    self.tbB.Size = Size(360, 25)
    self.tbB.Location = Point(100, 60)
    self.tbB.ValueChanged += self.updatetboxB
    self.Controls.Add(self.tbB)

    # Add a couch label
    labelC = Label()
    labelC.Text = 'Please select a couch angle in DEG [-90:+90].'
    labelC.Location = Point(15, 115)
    labelC.AutoSize = True
    self.Controls.Add(labelC)

    # Add a text box that to write the desired couch angle
    self.tboxC = TextBox()
    self.tboxC.Location = Point(15, 160)
    self.tboxC.Width = 55
    self.tboxC.Text = "0"
    self.tboxC.KeyDown += self.on_enter
    self.Controls.Add(self.tboxC)

    # Add a trackbar to slide to the desired couch angle
    self.tbC = TrackBar()
    self.tbC.TickStyle = TickStyle.Both
    self.tbC.TickFrequency = 5
    self.tbC.Minimum = -90
    self.tbC.Maximum = 90
    self.tbC.Value = 0
    self.tbC.Size = Size(360, 25)
    self.tbC.Location = Point(100, 160)
    self.tbC.ValueChanged += self.updatetboxC
    self.Controls.Add(self.tbC)

    # Add button to press Apply
    button = Button()
    button.Text = 'Apply'
    button.AutoSize = True
    button.Location = Point(15, 225)
    button.Click += self.apply_button_clicked
    self.Controls.Add(button)

    button2 = Button()
    button2.Text = 'Exit'
    button2.AutoSize = True
    button2.Location = Point(375, 225)
    button2.Click += self.exit_button_clicked
    self.Controls.Add(button2)

  def on_enter(self, sender, args):
    key = args.KeyCode
    if key == Keys.Enter:
      self.transform()

  def updatetboxB(self,sender,event):
    self.tboxB.Text = str(self.tbB.Value)
    self.transform()

  def updatetboxC(self, sender, event):
    self.tboxC.Text = str(self.tbC.Value)
    self.transform()

  def exit_button_clicked(self, sender, event):
    self.Close()

  def apply_button_clicked(self, sender, event):
    # Method invoked when the button is clicked
    self.transform()

  def transform(self):
    ba = self.tboxB.Text
    ca = self.tboxC.Text

    #Sanity check
    ok = True
    if ba=="" or float(ba)<self.tbB.Minimum:
      ba=str(int(self.tbB.Minimum))
      self.tboxB.Text = ba
      ok = False
    if ca=="" or float(ca)<self.tbC.Minimum:
      ca = str(int(self.tbC.Minimum))
      self.tboxC.Text = ca
      ok = False
    if float(ba) > self.tbB.Maximum:
      ba = str(int(self.tbB.Maximum))
      self.tboxB.Text = ba
      ok = False
    if float(ca) > self.tbC.Maximum:
      ca = str(int(self.tbC.Maximum))
      self.tboxC.Text = ca
      ok = False

    self.update_sliders()
    if ok==True:
      global gangle
      global cangle
      global oldgangle
      global oldcangle
      # Transform the models
      oldgangle = gangle
      oldcangle = cangle
      gangle = radians(float(ba))
      cangle = radians(float(ca))
      transform_models()

  def update_sliders(self):
    #Without emitting, to avoid inf loop
    self.tbB.ValueChanged -= self.updatetboxB
    self.tbB.Value = float(self.tboxB.Text)
    self.tbB.ValueChanged += self.updatetboxB
    self.tbC.ValueChanged -= self.updatetboxC
    self.tbC.Value = float(self.tboxC.Text)
    self.tbC.ValueChanged += self.updatetboxC

def tune_models():
  #await_user_input('Import finished. You can check your model now. If you want to change the beam angle, click on Resume again. If you want to exit and remove the ROIs, introduce a negative beam angle.')
  aform = SelectAngleForm()
  Application.Run(aform)
  remove_models()


def transform_models():
  for part in machine.parts:
    roiName = part.name
    if roiName=='Nozzle':
      b = -cs*(oldcangle+c0)
      b2 = cs*(cangle+c0)
      d = gs*(gangle - oldgangle) #g0 cancels
      case.PatientModel.RegionsOfInterest[roiName].TransformROI3D(Examination=examination, TransformationMatrix= {
        'M11':cos(d)*cos(b)*cos(b2)-sin(b)*sin(b2), 'M12':-sin(d)*cos(b2), 'M13':-cos(d)*sin(b)*cos(b2)-cos(b)*sin(b2), 'M14':iso.x-iso.x*(cos(d)*cos(b)*cos(b2)-sin(b)*sin(b2))+iso.y*sin(d)*cos(b2)+iso.z*(cos(d)*sin(b)*cos(b2)+cos(b)*sin(b2)),
        'M21':sin(d)*cos(b)                       , 'M22': cos(d)        , 'M23':-sin(d)*sin(b)                       , 'M24':iso.y-iso.x* sin(d)*cos(b)                        -iso.y*cos(d)        +iso.z* sin(d)*sin(b)                        ,
        'M31':cos(d)*cos(b)*sin(b2)+sin(b)*cos(b2), 'M32':-sin(d)*sin(b2), 'M33':-cos(d)*sin(b)*sin(b2)+cos(b)*cos(b2), 'M34':iso.z-iso.x*(cos(d)*cos(b)*sin(b2)+sin(b)*cos(b2))+iso.y*sin(d)*sin(b2)+iso.z*(cos(d)*sin(b)*sin(b2)-cos(b)*cos(b2)),
        'M41':0                                   , 'M42': 0             , 'M43': 0                                   , 'M44':1                                                                                                                   })
  #await_user_input('Transformation finished. If you want to change the beam angle, click on Resume Script again. If you want to remove the ROIs, introduce a negative beam angle.')


def remove_models():
  for part in machine.parts:
    # delete ROI
    roiName = part.name
    case.PatientModel.RegionsOfInterest[roiName].DeleteRoi()


#Initialization. Variables below are global
form = SelectMachineForm()
Application.Run(form)
mname = form.mach_name

case = get_current('Case')
examination = get_current('Examination')
machine = machines[mname]
structure_set = case.PatientModel.StructureSets[examination.Name]
orientation = examination.PatientPosition

# A rotation of the 3D model is needed to match the CT orientation depending on PatientPosition attribute
gantry_angle_offset = {'HFS':180, 'FFS':180,'HFP':0  , 'FFP': 0}
couch_angle_offset =  {'HFS':180, 'FFS':  0,'HFP':180, 'FFP': 0}
gantry_direction =    {'HFS': -1, 'FFS': -1,'HFP':-1 , 'FFP':-1}
couch_direction =     {'HFS': -1, 'FFS': -1,'HFP': 1 , 'FFP': 1}
g0 = radians(gantry_angle_offset[orientation])
c0 = radians(couch_angle_offset[orientation])
gs = gantry_direction[orientation]
cs = couch_direction[orientation]

poi_type = 'Isocenter'
poi_lst = [r.Type for r in case.PatientModel.PointsOfInterest]
#check if POI already defined, if not, wait until defined, then continue
while poi_type not in poi_lst:
  await_user_input('Please click OK and define an "'+poi_type+'" POI, then click on Play Script')
  poi_lst = [r.Type for r in case.PatientModel.PointsOfInterest]

iso = structure_set.PoiGeometries[poi_lst.index(poi_type)].Point

#Remove previous ROIs if already defined, e.g. if previous program instance crashed or script was stopped
roi_lst = [r.Name for r in case.PatientModel.RegionsOfInterest]
for part in machine.parts:
  # create ROI
  roiName = part.name
  if roiName in roi_lst:
    await_user_input('Confirm deletion of preexisting ROI "'+ roiName + '" by clicking on Resume Script. Otherwise click Stop Script.')
    case.PatientModel.RegionsOfInterest[roiName].DeleteRoi()

#Create first model at angles 0,0
gangle = 0
cangle = 0
oldgangle = 0
oldcangle = 0

for part in machine.parts:
  # create ROI
  roiName = part.name
  roiColor = part.color
  roiType = 'Support'
  fileName = machine.path+part.file
  case.PatientModel.CreateRoi(Name=roiName, Color=roiColor, Type=roiType)
  # import mesh from file
  geo = case.PatientModel.StructureSets[examination.Name].RoiGeometries[roiName]
  if part.name=='Nozzle':
    a = gs*(gangle+g0)
    b = cs*(cangle+c0)
    geo.ImportRoiGeometryFromSTL(FileName=fileName, UnitInFile='Millimeter',
    TransformationMatrix={'M11':cos(a)*cos(b), 'M12':-sin(a)*cos(b), 'M13':-sin(b), 'M14':iso.x,
                               'M21':sin(a)       , 'M22': cos(a)       , 'M23': 0     , 'M24':iso.y,
                               'M31':cos(a)*sin(b), 'M32':-sin(a)*sin(b), 'M33': cos(b), 'M34':iso.z,
                               'M41':0            , 'M42':0             , 'M43': 0     , 'M44':1          })

ScriptClient.AppUtil.RunInNewThread(tune_models())
