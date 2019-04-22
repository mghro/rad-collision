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

from math import cos, sin, radians, degrees, hypot, atan2
import sys
import re
import itertools
import ScriptClient
from connect import *
import clr
clr.AddReference("System.Windows.Forms")
clr.AddReference("System.Drawing")
from System.Windows.Forms import Application, Form, Label, ComboBox, Button, TextBox, TrackBar, FormStartPosition, TickStyle, Keys, CheckBox
from System.Drawing import Point, Size


class Part:

    def __init__(self, name, file, color, active, movex=True, movey=True, movez=True):
        self.name = name
        self.file = file
        self.color = color
        self.active = active
        self.moveX = movex
        self.moveY = movey
        self.moveZ = movez


class Machine:

    def __init__(self, name, path, parts):
        self.name = name
        self.path = path
        self.parts = parts


# agility = Machine("Agility", "\\\\Client\\H$\\STL parts\\",[Part("Nozzle","RotatingHeads.stl","Blue") ])
agility = Machine("Elekta Agility", "Y:\\STL parts\\Elekta Agility\\",
                  [Part("Nozzle", "RotatingHeads.stl", "Blue", True),
                   Part("Electron cone 15cm x 15cm", "ElectronApplicator.stl", "Blue", False),
                   Part("Flat panel left", "FlatPanelExtensionLeft.stl", "Blue", False),
                   Part("Flat panel bottom", "FlatPanelExtension.stl", "Blue", False),
                   Part("CBCT source", "CBCTsource.stl", "Blue", False)
                   ]
                  )
proteus = Machine("iba Proteus", "Y:\\STL parts\\iba Proteus\\", [Part("Nozzle", "NozzleWithSmallSnout.stl", "Blue", True)])
trilogy = Machine("Varian IX Trilogy", "Y:\\STL parts\\Varian IX Trilogy\\", [Part("Nozzle", "RotatingHeads.stl", "Blue", True)])
truebeam = Machine("Varian TrueBeam",   "Y:\\STL parts\\Varian TrueBeam\\", [Part("Nozzle", "RotatingHeads.stl", "Blue", True)])

linacs = {agility.name: agility, proteus.name: proteus, trilogy.name: trilogy, truebeam.name: truebeam}

evo = Machine("Hexapod Evo", "Y:\\STL parts\\Hexapod Evo\\",
              [Part("Couch", "Hexapod.stl", "Green", True),
               Part("Couch base top", "CouchBaseFwd.stl", "Green", False),
               Part("Couch base middle", "CouchBase.stl", "Green", False, True, True, False),
               Part("Couch base bottom", "ScissorsTop.stl", "Green", False, False, True, False),
               Part("Couch base support", "ScissorsBottom.stl", "Green", False, False, False, False),
               Part("Right control", "RightCouchControl.stl", "Green", False, True, True, False),
               Part("Left control", "LeftCouchControl.stl", "Green", False, True, True, False),
               Part("Extension", "CouchExtension.stl", "Green", False),
               Part("Head support", "HeadAdapter.stl", "Green", False),
               ]
              )

couches = {evo.name: evo}


# Define Form class that will prompt the user to select an element
# from a list via a combo box
class SelectListForm(Form):

    def __init__(self, list, description):
        # Set the size of the form
        self.name = None
        self.Size = Size(300, 200)
        # Set title of the form
        self.Text = 'Select '+description

        # Add a label
        label = Label()
        label.Text = 'Please select desired '+description
        label.Location = Point(15, 15)
        label.AutoSize = True
        self.Controls.Add(label)

        # Add a ComboBox that will display the Machines to select
        # Define the items to show up in the combobox
        self.combobox = ComboBox()
        self.combobox.DataSource = list.keys()
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

        self.AcceptButton = button

    def ok_button_clicked(self, _sender, _event):
        # Method invoked when the button is clicked
        # Save the selected machine name
        self.name = self.combobox.SelectedValue
        # Close the form
        self.Close()


# Define Form class that will prompt the user to select which
# parts to draw
class SelectPartsForm(Form):

    def __init__(self, group):
        # Set the size of the form
        self.Size = Size(300, 500)
        # Set title of the form
        self.Text = 'Select Parts'
        self.group = group

        # Add a label
        label = Label()
        label.Text = 'Please select parts to include'
        label.Location = Point(15, 15)
        label.AutoSize = True
        self.Controls.Add(label)

        # Add a CheckBox for each part to activate or deactivate
        i = 0
        for part in self.group.parts:
            part.cb = CheckBox()
            part.cb.Location = Point(15, 60+i*20)
            part.cb.Text = part.name
            part.cb.Width = 275
            part.cb.Checked = part.active
            self.Controls.Add(part.cb)
            i += 1

        # Add button to press OK and close the form
        button = Button()
        button.Text = 'OK'
        button.AutoSize = True
        button.Location = Point(15, 400)
        button.Click += self.ok_button_clicked
        self.Controls.Add(button)

        self.AcceptButton = button

    def ok_button_clicked(self, _sender, _event):
        # Method invoked when the button is clicked
        # Save the selected machine name
        for part in self.group.parts:
            part.active = part.cb.Checked
            part.cb = None
        # Close the form
        self.Close()


class SelectAngleForm(Form):
    def __init__(self):
        # Set the size of the form
        self.StartupPosition = FormStartPosition.Manual
        self.Location = Point(500, 15)
        self.Size = Size(500, 475)
        # Set title of the form
        self.Text = 'Tune 3D model positions'
        self.TopMost = True

        # Add a beam label
        label_b = Label()
        label_b.Text = 'Please select a beam angle in DEG [0:360].'
        label_b.Location = Point(15, 15)
        label_b.AutoSize = True
        self.Controls.Add(label_b)

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
        self.tbB.TickFrequency = 10
        self.tbB.Minimum = 0
        self.tbB.Maximum = 360
        self.tbB.Value = 0
        self.tbB.Size = Size(360, 25)
        self.tbB.Location = Point(100, 60)
        self.tbB.ValueChanged += self.updatetbox_b
        self.Controls.Add(self.tbB)

        # Add a couch angle label
        label_c = Label()
        label_c.Text = 'Please select a couch angle in DEG [-90:+90].'
        label_c.Location = Point(15, 115)
        label_c.AutoSize = True
        self.Controls.Add(label_c)

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
        self.tbC.ValueChanged += self.updatetbox_c
        self.Controls.Add(self.tbC)

        # Add a couch xyz label
        label_xyz = Label()
        label_xyz.Text = 'Please select XYZ couch position  [mm].'
        label_xyz.Location = Point(15, 215)
        label_xyz.AutoSize = True
        self.Controls.Add(label_xyz)

        # Add a text box to write the desired couch x position
        self.tboxX = TextBox()
        self.tboxX.Location = Point(15, 260)
        self.tboxX.Width = 55
        self.tboxX.Text = "0"
        self.tboxX.KeyDown += self.on_enter
        self.Controls.Add(self.tboxX)

        # Add a trackbar to slide to the desired couch x position
        self.tbX = TrackBar()
        self.tbX.TickStyle = TickStyle.Both
        self.tbX.TickFrequency = 10
        self.tbX.Minimum = -100
        self.tbX.Maximum = 100
        self.tbX.Value = 0
        self.tbX.Size = Size(360, 25)
        self.tbX.Location = Point(100, 260)
        self.tbX.ValueChanged += self.updatetbox_x
        self.Controls.Add(self.tbX)

        # Add a text box to write the desired couch y position
        self.tboxY = TextBox()
        self.tboxY.Location = Point(15, 300)
        self.tboxY.Width = 55
        self.tboxY.Text = "0"
        self.tboxY.KeyDown += self.on_enter
        self.Controls.Add(self.tboxY)

        # Add a trackbar to slide to the desired couch y position
        self.tbY = TrackBar()
        self.tbY.TickStyle = TickStyle.Both
        self.tbY.TickFrequency = 50
        self.tbY.Minimum = -250
        self.tbY.Maximum = 250
        self.tbY.Value = 0
        self.tbY.Size = Size(360, 25)
        self.tbY.Location = Point(100, 300)
        self.tbY.ValueChanged += self.updatetbox_y
        self.Controls.Add(self.tbY)

        # Add a text box to write the desired couch z position
        self.tboxZ = TextBox()
        self.tboxZ.Location = Point(15, 340)
        self.tboxZ.Width = 55
        self.tboxZ.Text = "0"
        self.tboxZ.KeyDown += self.on_enter
        self.Controls.Add(self.tboxZ)

        # Add a trackbar to slide to the desired couch y position
        self.tbZ = TrackBar()
        self.tbZ.TickStyle = TickStyle.Both
        self.tbZ.TickFrequency = 100
        self.tbZ.Minimum = -500
        self.tbZ.Maximum = 500
        self.tbZ.Value = 0
        self.tbZ.Size = Size(360, 25)
        self.tbZ.Location = Point(100, 340)
        self.tbZ.ValueChanged += self.updatetbox_z
        self.Controls.Add(self.tbZ)

        # Add button to press Apply
        button = Button()
        button.Text = 'Apply'
        button.AutoSize = True
        button.Location = Point(15, 390)
        button.Click += self.apply_button_clicked
        self.Controls.Add(button)

        button2 = Button()
        button2.Text = 'Exit'
        button2.AutoSize = True
        button2.Location = Point(375, 390)
        button2.Click += self.exit_button_clicked
        self.Controls.Add(button2)

    def on_enter(self, _sender, args):
        key = args.KeyCode
        if key == Keys.Enter:
            self.transform()

    def updatetbox_b(self, _sender, _event):
        self.tboxB.Text = str(self.tbB.Value)
        self.transform()

    def updatetbox_c(self, _sender, _event):
        self.tboxC.Text = str(self.tbC.Value)
        self.transform()

    def updatetbox_x(self, _sender, _event):
        self.tboxX.Text = str(self.tbX.Value)
        self.transform()

    def updatetbox_y(self, _sender, _event):
        self.tboxY.Text = str(self.tbY.Value)
        self.transform()

    def updatetbox_z(self, _sender, _event):
        self.tboxZ.Text = str(self.tbZ.Value)
        self.transform()

    def exit_button_clicked(self, _sender, _event):
        self.Close()

    def apply_button_clicked(self, _sender, _event):
        # Method invoked when the button is clicked
        self.transform()

    def transform(self):
        ba = self.tboxB.Text
        ca = self.tboxC.Text
        x = self.tboxX.Text
        y = self.tboxY.Text
        z = self.tboxZ.Text

        # Sanity check
        ok = True
        if ba == "" or float(ba) < self.tbB.Minimum:
            ba = str(int(self.tbB.Minimum))
            self.tboxB.Text = ba
            ok = False
        if ca == "" or float(ca) < self.tbC.Minimum:
            ca = str(int(self.tbC.Minimum))
            self.tboxC.Text = ca
            ok = False
        if x == "" or float(x) < self.tbX.Minimum:
            x = str(int(self.tbX.Minimum))
            self.tboxX.Text = x
            ok = False
        if y == "" or float(y) < self.tbY.Minimum:
            y = str(int(self.tbY.Minimum))
            self.tboxY.Text = y
            ok = False
        if z == "" or float(z) < self.tbZ.Minimum:
            y = str(int(self.tbZ.Minimum))
            self.tboxZ.Text = z
            ok = False
        if float(ba) > self.tbB.Maximum:
            ba = str(int(self.tbB.Maximum))
            self.tboxB.Text = ba
            ok = False
        if float(ca) > self.tbC.Maximum:
            ca = str(int(self.tbC.Maximum))
            self.tboxC.Text = ca
            ok = False
        if float(x) > self.tbX.Maximum:
            x = str(int(self.tbX.Maximum))
            self.tboxX.Text = x
            ok = False
        if float(y) > self.tbY.Maximum:
            y = str(int(self.tbY.Maximum))
            self.tboxY.Text = y
            ok = False
        if float(z) > self.tbZ.Maximum:
            z = str(int(self.tbZ.Maximum))
            self.tboxZ.Text = z
            ok = False

        self.update_sliders()
        if ok:
            global gangle
            global cangle
            global oldgangle
            global oldcangle
            global cx
            global oldcx
            global cy
            global oldcy
            global cz
            global oldcz
            # Transform the models
            oldgangle = gangle
            oldcangle = cangle
            oldcx = cx
            oldcy = cy
            oldcz = cz
            gangle = radians(float(ba))
            cangle = radians(float(ca))
            cx = float(x)
            cy = float(y)
            cz = float(z)
            #Convert to cm
            cx /= 10.
            cy /= 10.
            cz /= 10.
            transform_models()

    def update_sliders(self):
        # Without emitting, to avoid inf loop
        self.tbB.ValueChanged -= self.updatetbox_b
        self.tbB.Value = float(self.tboxB.Text)
        self.tbB.ValueChanged += self.updatetbox_b
        self.tbC.ValueChanged -= self.updatetbox_c
        self.tbC.Value = float(self.tboxC.Text)
        self.tbC.ValueChanged += self.updatetbox_c
        self.tbX.ValueChanged -= self.updatetbox_x
        self.tbX.Value = float(self.tboxX.Text)
        self.tbX.ValueChanged += self.updatetbox_x
        self.tbY.ValueChanged -= self.updatetbox_y
        self.tbY.Value = float(self.tboxY.Text)
        self.tbY.ValueChanged += self.updatetbox_y
        self.tbZ.ValueChanged -= self.updatetbox_z
        self.tbZ.Value = float(self.tboxZ.Text)
        self.tbZ.ValueChanged += self.updatetbox_z


def tune_models():
    # await_user_input('Import finished. You can check your model now. If you want to change the beam angle,
    # click on Resume again. If you want to exit and remove the ROIs, introduce a negative beam angle.')
    aform = SelectAngleForm()
    Application.Run(aform)
    remove_models()


def transform_models():
    for part in machine.parts:
        if part.active:
            if abs(cangle-oldcangle) > 0 or abs(gangle-oldgangle) > 0:
                roi_name = part.name
                b = -cs*(oldcangle+c0)
                b2 = cs*(cangle+c0)
                d = gs*(gangle - oldgangle)  # g0 cancels
                case.PatientModel.RegionsOfInterest[roi_name].TransformROI3D(Examination=examination, TransformationMatrix={
                    'M11': cos(d)*cos(b)*cos(b2)-sin(b)*sin(b2), 'M12': -sin(d)*cos(b2), 'M13': -cos(d)*sin(b)*cos(b2)-cos(b)*sin(b2), 'M14': iso.x-iso.x*(cos(d)*cos(b)*cos(b2)-sin(b)*sin(b2))+iso.y*sin(d)*cos(b2)+iso.z*(cos(d)*sin(b)*cos(b2)+cos(b)*sin(b2)),
                    'M21': sin(d)*cos(b)                       , 'M22':  cos(d)        , 'M23': -sin(d)*sin(b)                       , 'M24': iso.y-iso.x* sin(d)*cos(b)                        -iso.y*cos(d)        +iso.z* sin(d)*sin(b)                        ,
                    'M31': cos(d)*cos(b)*sin(b2)+sin(b)*cos(b2), 'M32': -sin(d)*sin(b2), 'M33': -cos(d)*sin(b)*sin(b2)+cos(b)*cos(b2), 'M34': iso.z-iso.x*(cos(d)*cos(b)*sin(b2)+sin(b)*cos(b2))+iso.y*sin(d)*sin(b2)+iso.z*(cos(d)*sin(b)*sin(b2)-cos(b)*cos(b2)),
                    'M41': 0                                   , 'M42':  0             , 'M43':  0                                   , 'M44': 1                                                                                                                   })
    # await_user_input('Transformation finished. If you want to change the beam angle, click on Resume Script again. If you want to remove the ROIs, introduce a negative beam angle.')
    for part in couch.parts:
        if part.active:
            roi_name = part.name
            dx = cx-oldcx
            dy = cy-oldcy
            dz = cz-oldcz
            if not part.moveX:
                dx = 0
            if not part.moveY:
                dy = 0
            if not part.moveZ:
                dz = 0
            if abs(dx) > 0 or abs(dy) > 0 or abs(dz) > 0:
                case.PatientModel.RegionsOfInterest[roi_name].TransformROI3D(Examination=examination, TransformationMatrix={
                    'M11': 1, 'M12': 0, 'M13': 0, 'M14': dx,
                    'M21': 0, 'M22': 1, 'M23': 0, 'M24': dy,
                    'M31': 0, 'M32': 0, 'M33': 1, 'M34': dz,
                    'M41': 0, 'M42': 0, 'M43': 0, 'M44': 1})


def remove_models():
    for part in itertools.chain(machine.parts, couch.parts):
        if part.active:
            # delete ROI
            roi_name = part.name
            case.PatientModel.RegionsOfInterest[roi_name].DeleteRoi()


def main():
    # Initialization. Variables below are global
    form = SelectListForm(linacs, "LINAC")
    Application.Run(form)
    global machine
    machine = linacs[form.name]

    pform = SelectPartsForm(machine)
    Application.Run(pform)

    cform = SelectListForm(couches, "patient couch")
    Application.Run(cform)
    global couch
    couch = couches[cform.name]

    pcform = SelectPartsForm(couch)
    Application.Run(pcform)

    global case
    case = get_current('Case')
    global examination
    examination = get_current('Examination')
    structure_set = case.PatientModel.StructureSets[examination.Name]
    orientation = examination.PatientPosition

    # A rotation of the 3D model is needed to match the CT orientation depending on PatientPosition attribute
    gantry_angle_offset = {'HFS': 180, 'FFS': 180, 'HFP':   0, 'FFP':  0}
    couch_angle_offset =  {'HFS': 180, 'FFS':   0, 'HFP': 180, 'FFP':  0}
    gantry_direction =    {'HFS':  -1, 'FFS':  -1, 'HFP':  -1, 'FFP': -1}
    couch_direction =     {'HFS':  -1, 'FFS':  -1, 'HFP':   1, 'FFP':  1}
    global g0, c0, gs, cs
    g0 = radians(gantry_angle_offset[orientation])
    c0 = radians(couch_angle_offset[orientation])
    gs = gantry_direction[orientation]
    cs = couch_direction[orientation]

    poi_type = 'Isocenter'
    poi_lst = [r.Type for r in case.PatientModel.PointsOfInterest]
    # check if POI already defined, if not, wait until defined, then continue
    while poi_type not in poi_lst:
        await_user_input('Please click OK and define an "'+poi_type+'" POI, then click on Play Script')
        poi_lst = [r.Type for r in case.PatientModel.PointsOfInterest]

    global iso
    iso = structure_set.PoiGeometries[poi_lst.index(poi_type)].Point

    # Create first model at angles 0,0.  These below are global variables.
    global gangle, cangle, oldgangle, oldcangle
    global cx, cy, cz, oldcx, oldcy, oldcz
    gangle = 0
    cangle = 0
    oldgangle = 0
    oldcangle = 0
    cx = 0
    cy = 0
    cz = 0
    oldcx = 0
    oldcy = 0
    oldcz = 0

    # Remove previous ROIs if already defined, e.g. if previous program instance crashed or script was stopped
    roi_lst = [r.Name for r in case.PatientModel.RegionsOfInterest]
    for part in itertools.chain(machine.parts, couch.parts):
        if part.active:
            # create ROI
            roi_name = part.name
            if roi_name in roi_lst:
                await_user_input('Confirm deletion of preexisting ROI "' + roi_name + '" by clicking on Resume Script. Otherwise click Stop Script.')
                case.PatientModel.RegionsOfInterest[roi_name].DeleteRoi()

    for part in machine.parts:
        if part.active:
            # create ROI
            roi_name = part.name
            roi_color = part.color
            roi_type = 'Support'
            file_name = machine.path+part.file
            case.PatientModel.CreateRoi(Name=roi_name, Color=roi_color, Type=roi_type)
            # import mesh from file
            geo = structure_set.RoiGeometries[roi_name]
            a = gs*(gangle+g0)
            b = cs*(cangle+c0)
            geo.ImportRoiGeometryFromSTL(FileName=file_name, UnitInFile='Millimeter',
                                         TransformationMatrix={'M11': cos(a)*cos(b), 'M12': -sin(a)*cos(b), 'M13': -sin(b), 'M14': iso.x,
                                                               'M21': sin(a)       , 'M22':  cos(a)       , 'M23':       0, 'M24': iso.y,
                                                               'M31': cos(a)*sin(b), 'M32': -sin(a)*sin(b), 'M33':  cos(b), 'M34': iso.z,
                                                               'M41':             0, 'M42':              0, 'M43':       0, 'M44':     1})

    # Get list of couches defined before the script
    couch_lst = [r.Name for r in case.PatientModel.RegionsOfInterest if r.Type == 'Support' if re.search('couch', r.Name, re.IGNORECASE)]

    for part in couch.parts:
        if part.active:
            # create ROI
            roi_name = part.name
            roi_color = part.color
            roi_type = 'Support'
            file_name = couch.path+part.file
            case.PatientModel.CreateRoi(Name=roi_name, Color=roi_color, Type=roi_type)
            # import mesh from file
            geo = structure_set.RoiGeometries[roi_name]
            a = gs*g0
            b = cs*c0
            geo.ImportRoiGeometryFromSTL(FileName=file_name, UnitInFile='Millimeter',
                                         TransformationMatrix={'M11': cos(a)*cos(b), 'M12': -sin(a)*cos(b), 'M13': -sin(b), 'M14': iso.x,
                                                               'M21': sin(a)       , 'M22':  cos(a)       , 'M23':       0, 'M24': iso.y,
                                                               'M31': cos(a)*sin(b), 'M32': -sin(a)*sin(b), 'M33':  cos(b), 'M34': iso.z,
                                                               'M41':             0, 'M42':              0, 'M43':       0, 'M44':     1})

    # If there is a Couch model and couch contour, recenter couch parts
    couch_models = [c.name for c in couch.parts if c.active if re.search('couch', c.name, re.IGNORECASE)]
    if len(couch_models) > 0:
        model_name = couch_models[0]
        geom = structure_set.RoiGeometries[model_name]
        if len(couch_lst) > 0:
            roi_name = couch_lst[0]
            geo = structure_set.RoiGeometries[roi_name]

            rb = geo.GetBoundingBox()
            rx = (rb[1].x + rb[0].x) / 2
            ry = (rb[1].y + rb[0].y) / 2
            rz = rb[0].z

            mb = geom.GetBoundingBox()
            mx = (mb[1].x + mb[0].x) / 2
            my = (mb[1].y + mb[0].y) / 2
            mz = mb[0].z

            dx0 = rx-mx
            dy0 = ry-my
            dz0 = rz-mz
            for part in couch.parts:
                if part.active:
                    roi_name = part.name
                    dx = dx0
                    dy = dy0
                    dz = dz0
                    if not part.moveX:
                        dx = 0
                    if not part.moveY:
                        dy = 0
                    if not part.moveZ:
                        dz = 0
                    case.PatientModel.RegionsOfInterest[roi_name].TransformROI3D(Examination=examination, TransformationMatrix={
                        'M11': 1, 'M12': 0, 'M13': 0, 'M14': dx,
                        'M21': 0, 'M22': 1, 'M23': 0, 'M24': dy,
                        'M31': 0, 'M32': 0, 'M33': 1, 'M34': dz,
                        'M41': 0, 'M42': 0, 'M43': 0, 'M44': 1})

    ScriptClient.AppUtil.RunInNewThread(tune_models())


if __name__ == '__main__':
    main()
