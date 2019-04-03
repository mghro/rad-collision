# Fernando Hueso-Gonzalez
# Massachusetts General Hospital and Harvard Medical School
#
# Illustrates how you can load a 3D model of the LINAC into RayStation for collision detection

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

agility = Machine("Agility", "\\\\Client\\D$\\Dropbox (Partners Healthcare)\\LINAC\\Elekta Agility\\STL parts\\",[Part("Nozzle","RotatingHeads.stl","Blue") ])

machines = {agility.name : agility}

import math
import sys

from connect import *

import clr
clr.AddReference("System.Windows.Forms")
clr.AddReference("System.Drawing")

from System.Windows.Forms import Application, Form, Label, ComboBox, Button, TextBox, TrackBar
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
    self.Size = Size(500, 200)
    # Set title of the form
    self.Text = 'Select Beam Angle'

    # Add a label
    label = Label()
    label.Text = 'Please select a beam angle in DEG [0-360]. Type -1 to exit.'
    label.Location = Point(15, 15)
    label.AutoSize = True
    self.Controls.Add(label)

    # Add a text box that to write the desired angle
    self.tbox = TextBox()
    self.tbox.Location = Point(15, 60)
    self.tbox.Width = 55
    self.Controls.Add(self.tbox)

    self.tb = TrackBar()
    #self.tb.TickStyle = TickStyle.Both
    self.tb.Minimum = 0
    self.tb.Maximum = 360
    self.tb.Size = Size(360, 25)
    self.tb.Location = Point(100, 60)
    self.tb.ValueChanged += self.updatetbox
    self.Controls.Add(self.tb)

    # Add button to press OK and close the form
    button = Button()
    button.Text = 'OK'
    button.AutoSize = True
    button.Location = Point(15, 100)
    button.Click += self.ok_button_clicked
    self.Controls.Add(button)

  def updatetbox(self,sender,event):
    self.tbox.Text = str(self.tb.Value)

  def ok_button_clicked(self, sender, event):
    # Method invoked when the button is clicked
    # Save the selected machine name
    self.beam_angle = self.tbox.Text
    # Close the form
    self.Close()

form = SelectMachineForm()
Application.Run(form)
mname = form.mach_name

case = get_current('Case')
examination = get_current('Examination')
machine = machines[mname]
structure_set = case.PatientModel.StructureSets[examination.Name]

poi_name = 'Iso'
poi_lst = [r.Name for r in case.PatientModel.PointsOfInterest]
#check if POI already defined, if not, wait until defined, then continue
while poi_name not in poi_lst:
  await_user_input('Please click OK and define an "'+poi_name+'" POI, then click on Play Script')
  poi_lst = [r.Name for r in case.PatientModel.PointsOfInterest]

isopoint = structure_set.PoiGeometries[poi_name].Point

aform = SelectAngleForm()
Application.Run(aform)
gangle = math.radians(float(aform.beam_angle))
if gangle < 0:
  sys.exit()

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
    geo.ImportRoiGeometryFromSTL(FileName=fileName, UnitInFile='Millimeter',
        TransformationMatrix={'M11':math.cos(-gangle), 'M12':-math.sin(-gangle), 'M13':0, 'M14':isopoint.x,
                              'M21':math.sin(-gangle), 'M22': math.cos(-gangle), 'M23':0, 'M24':isopoint.y,
                              'M31':0                , 'M32':0                 , 'M33':1, 'M34':isopoint.z,
                              'M41':0                , 'M42':0                 , 'M43':0, 'M44':1          })

await_user_input('Import finished. You can check your model now. If you want to change the beam angle, click on Resume again. If you want to exit, introduce a negative beam angle. The LINAC will be removed from the ROI list.')

while gangle >= 0:
  aform = SelectAngleForm()
  Application.Run(aform)
  oldgangle = gangle
  gangle = math.radians(float(aform.beam_angle))
  if gangle>=0:
    for part in machine.parts:
      if part.name=='Nozzle':
        tangle = gangle - oldgangle
        x = isopoint.x
        y = isopoint.y
        h = math.hypot(x,y)
        alpha = math.atan2(y,x)
        case.PatientModel.RegionsOfInterest[part.name].TransformROI3D(Examination=examination,
          TransformationMatrix={'M11':math.cos(-tangle), 'M12':-math.sin(-tangle), 'M13':0, 'M14':x-h*math.cos(alpha-tangle),
                                'M21':math.sin(-tangle), 'M22': math.cos(-tangle), 'M23':0, 'M24':y-h*math.sin(alpha-tangle),
                                'M31':0               , 'M32':0                , 'M33':1, 'M34':0,
                                'M41':0               , 'M42':0                , 'M43':0, 'M44':1          })
    #https://www.tutorialspoint.com/computer_graphics/3d_transformation.htm
    await_user_input('Transformation finished. If you want to change the beam angle, click on Resume Script again. If you want to exit, introduce a negative beam angle. The LINAC will be removed from the ROI list.')

for part in machine.parts:
  # delete ROI
  roiName = part.name
  case.PatientModel.RegionsOfInterest[roiName].DeleteRoi()
