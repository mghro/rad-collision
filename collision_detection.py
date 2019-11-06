#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# collision_detection.py
#
# Fernando Hueso-González - fhuesogonzalez@mgh.harvard.edu
# Massachusetts General Hospital and Harvard Medical School
#
# Loads a 3D model of a radiotherapy treatment head and patient couch into RayStation for collision detection
#
# The 3D models are not part of the script, you need to ask your vendor to provide them as STL files, eg. under an NDA.
# They have to be stored in a folder visible by the RayStation server. The different STL parts of the treatment head
# must be stored within the same folder. The couch has to be stored in another folder with its corresponding subparts.
# For example:
# F:\\STL models\\
#                 LINAC\\
#                        Nozzle.stl
#                        XrayPanel.stl
#                 Couch\\
#                        Hexapod.stl
#                        HeadSupport.stl
# If the STL files are not stored on the server but on the client side, you have to specify the path as
# \\\\Client\\F\\STL models\\
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
# However, the 2D viewer looks different for HFS orientation or any other, as they maintain the CT acquisition orientation
# and move instead the labels R, L, A, P, S, I.
# Note also that the coordinates you see in the RayStation Viewer are not xyz, but RL, IS, PA, with RL = x, IS = z, PA = -y
#
# Let us define some rotation angles:
# a: rotation around z axis (from x to y axis), represented by a classical rotation matrix R_z, with R_z[row=1,column=2] = -sin(a)
# b: rotation around y axis (from x to z axis), represented by a classical rotation matrix R_y, with R_y[row=1,column=3] = -sin(b)
# g: IEC DICOM gantry rotation. When couch angle is zero, in the case of FFS, it is from y to x axis (opposite sign than a).
# c: IEC Couch support rotation. When gantry angle is zero, in the case of FFS, it is from x to z axis (same sign than b).
#
# It should be noted that the sign of the rotation angle depends on the patient orientation (FFS, etc.).
# For each patient orientation, we define the sign (gs,cs) between gantry and couch rotation (g,c) and the respective rotation
# about the DICOM patient axes (a,b). Also, the gantry and couch offset (g0,c0) are defined to rotate the 3D model in order
# to match the DICOM patient axes at the particular patient orientation.
#
# Coordinates of the 3D model
# It is a requisite that the origin of the STL 3D models is exactly at room isocenter.
# Also, the orientation should be such that, when opening the STL file with Meshlab, you are looking towards the treatment head
# as if you were standing in the treatment room in front of it (no couch support rotation), and the gantry is at zero degrees.
# You might need to cleanup your model to avoid an import error in RayStation. You can do so by opening it in Meshlab, then
# click on Unify duplicated Vertices, Ok, File, Export Mesh As, Select first extension on dropdown menu, and then rename as .stl
#
# The initial affine transformation will place the treatment head at gantry and couch angle g=0, c=0.
# Internally, we rotate the treatment head according to gantry angle offset g0,
# then simulate couch angle offset c0 via a negative couch rotation of the model,
# and finally the translation of the model to the isocenter (iso.x,iso.y,iso.z) in the CT patient
# resulting in TransformationMatrix M(iso.x,iso.y,iso.z,c,g) = T(iso.x,iso.y,iso.z) * R_y ( b = -(c+c0) ) * R_z ( a = -(g+g0) ) =
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
# M(iso.x,iso.y,iso.z,c2,g2) * R_z ( a = g+g0 ) * R_y ( b = c+c0 ) * T(-iso.x,-iso.y,-iso.z) =
# T(iso.x,iso.y,iso.z) * R_y ( b2 = -(c0+c2) ) * R_z ( a2 = -(g0+g2))  * R_z ( a = (g0+g) ) * R_y ( b = (c0+c) ) * T(-iso.x,-iso.y,-iso.z) =
# T(iso.x,iso.y,iso.z) * R_y ( b2 = -(c0+c2) ) * R_z ( d = g-g2) * R_y ( b = c ) * T(-iso.x,-iso.y,-iso.z) =
# resulting in TransformationMatrix D(iso.x,iso.y,iso.z,c2,g2,c,g) =
# {'M11':cos(d)*cos(b)*cos(b2)-sin(b)*sin(b2), 'M12':-sin(d)*cos(b2), 'M13':-cos(d)*sin(b)*cos(b2)-cos(b)*sin(b2), 'M14':iso.x-iso.x*(cos(d)*cos(b)*cos(b2)-sin(b)*sin(b2))+iso.y*sin(d)*cos(b2)+iso.z*(cos(d)*sin(b)*cos(b2)+cos(b)*sin(b2)),
#  'M21':sin(d)*cos(b)                       , 'M22': cos(d)        , 'M23':-sin(d)*sin(b)                       , 'M24':iso.y-iso.x* sin(d)*cos(b)                        -iso.y*cos(d)        +iso.z* sin(d)*sin(b)                        ,
#  'M31':cos(d)*cos(b)*sin(b2)+sin(b)*cos(b2), 'M32':-sin(d)*sin(b2), 'M33':-cos(d)*sin(b)*sin(b2)+cos(b)*cos(b2), 'M34':iso.z-iso.x*(cos(d)*cos(b)*sin(b2)+sin(b)*cos(b2))+iso.y*sin(d)*sin(b2)+iso.z*(cos(d)*sin(b)*sin(b2)-cos(b)*cos(b2)),
#  'M41':0                                   , 'M42': 0             , 'M43': 0                                   , 'M44':1                                                             }

# Import basic modules
from math import cos, sin, radians, degrees, sqrt, acos, atan2
import re
import itertools

# Import RayStation modules and WinForms for GUI
from connect import *
import clr
clr.AddReference("System.Windows.Forms")
clr.AddReference("System.Drawing")
from System.Windows.Forms import Application, Form, Label, ComboBox, Button, TextBox, TrackBar, FormStartPosition, TickStyle, Keys, CheckBox
from System.Drawing import Point, Size
from System.Threading import ThreadStart, Thread


class Part:
    """
    Class describing a 3D model file, that might be a part of the whole machine (treatment head or couch)
    """

    def __init__(self, name, filename, color, active, movex=True, movey=True, movez=True, scissor=False, retractable=False):
        """
        Initialization of the object
        :param name: the identifier name of the part, it must be unique as it is used as a key, e.g. Nozzle.
        :param filename: the name of the STL file within the folder you stored it
        :param color: the color of the ROI once the model is imported into Raystation
        :param active: flag to activate or deactivate the import of this specific part of the whole machine
        :param movex: flag to activate or deactivate the translation of this part along the x coordinate
        :param movey: flag to activate or deactivate the translation of this part along the y coordinate
        :param movez: flag to activate or deactivate the translation of this part along the z coordinate
        :param scissor: if this part is a scissor Robot
        :param retractable: flag to signal if this part is retractable, i.e. a snout or range shifter in the nozzle
        """
        self.name = name
        self.filename = filename
        self.color = color
        self.active = active
        self.moveX = movex
        self.moveY = movey
        self.moveZ = movez
        self.scissor = scissor
        self.retractable = retractable


class Machine:
    """
    Class grouping different Parts into the same machine, e.g. a treatment head or a couch
    """

    def __init__(self, name, path, parts):
        """
        Initialization of the object
        :param name: the identifier name of the part, it must be unique as it is used as a key, e.g. Elekta Agility.
        :param path: the path where all STL models of this machine are stored, namely the folder containing all subfiles (STL parts)
        :param parts: the array of Parts corresponding to this machine
        """
        self.name = name
        self.path = path
        self.parts = parts


class SelectListForm(Form):
    """
    Define Form generic class that will prompt the user to select an element
    from a list via a combo box
    """

    def __init__(self, lst, description):
        """
        :param: self the reference to the Form
        :param: lst the list of elements with keys
        :param: description the readable title of the list or category of elements
        """
        self.name = None
        self.Size = Size(300, 200)   # Set the size of the form
        self.Text = 'Select {}'.format(description)  # Set title of the form

        # Add a label
        label = Label()
        label.Text = 'Please select desired {}'.format(description)
        label.Location = Point(15, 15)
        label.AutoSize = True
        self.Controls.Add(label)

        # Add a ComboBox that will display the items of this list
        self.combobox = ComboBox()
        self.combobox.DataSource = lst.keys()
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

        self.AcceptButton = button  # Enter Key presses the OK button

    def ok_button_clicked(self, _sender, _event):
        """
        Method invoked when the OK button is clicked
        :param: self the reference to the Form
        :param: _sender  ignore
        :param: _event ignore
        """
        self.name = self.combobox.SelectedValue  # Save the selected element name
        self.Close()  # Close the form


class SelectPartsForm(Form):
    """
    Define Form class that will prompt the user to select which parts of the Machine to draw
    """

    def __init__(self, machine):
        """
        Form initialization
        :param: self reference to the Form
        :param: machine array of parts corresponding to this machine
        """

        self.Size = Size(300, 500)  # Set the size of the form
        self.Text = 'Select Parts'  # Set title of the form
        self.machine = machine

        # Add a label
        label = Label()
        label.Text = 'Please select parts to include'
        label.Location = Point(15, 15)
        label.AutoSize = True
        self.Controls.Add(label)

        # Add a CheckBox for each part to activate or deactivate each part separately
        for i, part in enumerate(self.machine.parts):
            part.cb = CheckBox()
            part.cb.Location = Point(15, 60+i*20)
            part.cb.Text = part.name
            part.cb.Width = 275
            part.cb.Checked = part.active
            self.Controls.Add(part.cb)

        # Add button to press OK and close the form
        button = Button()
        button.Text = 'OK'
        button.AutoSize = True
        button.Location = Point(15, 400)
        button.Click += self.ok_button_clicked
        self.Controls.Add(button)

        self.AcceptButton = button  # Pressing enter works like a click on the OK button

    def ok_button_clicked(self, _sender, _event):
        """
        Method invoked when the OK button is clicked. Parts active flag is updated
        :param: self the reference to the Form
        :param: _sender  ignore
        :param: _event ignore
        """
        for part in self.machine.parts:
            part.active = part.cb.Checked
            part.cb = None
        self.Close()  # Close the form


class TuneModelsForm(Form):
    """
    Main GUI form to move and rotate 3D models interactively
    """

    def __init__(self):
        """
        Form initialization
        :param: self reference to the Form
        """
        self.StartupPosition = FormStartPosition.Manual
        self.Location = Point(500, 15)
        self.Size = Size(500, 575 if extraction else 475)  # Set the size of the form
        self.Text = 'Tune 3D model positions'  # Set title of the form
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
        self.tbB.Location = Point(100, 60-5)
        self.tbB.ValueChanged += self.updatetbox_b
        self.Controls.Add(self.tbB)

        # Add a couch angle label
        label_c = Label()
        label_c.Text = 'Please select a couch angle in DEG [-90:+90].'
        label_c.Location = Point(15, 115)
        label_c.AutoSize = True
        self.Controls.Add(label_c)

        # Add a text box to write the desired couch angle
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
        self.tbC.Location = Point(100, 160-5)
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
        self.tbX.Location = Point(100, 260-5)
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
        self.tbY.Location = Point(100, 300-5)
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
        self.tbZ.Location = Point(100, 340-5)
        self.tbZ.ValueChanged += self.updatetbox_z
        self.Controls.Add(self.tbZ)

        lastpos = 340

        if extraction:
            # Add a retraction label
            label_ext = Label()
            label_ext.Text = 'Please select snout extraction  [mm].'
            label_ext.Location = Point(15, 395)
            label_ext.AutoSize = True
            self.Controls.Add(label_ext)

            # Add a text box to write the desired couch x position
            self.tboxE = TextBox()
            self.tboxE.Location = Point(15, 440)
            self.tboxE.Width = 55
            self.tboxE.Text = "0"
            self.tboxE.KeyDown += self.on_enter
            self.Controls.Add(self.tboxE)

            # Add a trackbar to slide to the desired couch x position
            self.tbE = TrackBar()
            self.tbE.TickStyle = TickStyle.Both
            self.tbE.TickFrequency = 40
            self.tbE.Minimum = 0
            self.tbE.Maximum = 800
            self.tbE.Value = 0
            self.tbE.Size = Size(360, 25)
            self.tbE.Location = Point(100, 440-5)
            self.tbE.ValueChanged += self.updatetbox_e
            self.Controls.Add(self.tbE)

            lastpos = 440

        # Add button to press Apply
        button = Button()
        button.Text = 'Apply'
        button.AutoSize = True
        button.Location = Point(15, lastpos+50)
        button.Click += self.apply_button_clicked
        self.Controls.Add(button)

        # Add button to press Flip in case of robot scissors
        if len(lsci) >= 2:
            button3 = Button()
            button3.Text = 'Flip'
            button3.AutoSize = True
            button3.Location = Point(125, lastpos+50)
            button3.Click += self.flip_button_clicked
            self.Controls.Add(button3)

        # Add button to press Exit
        button2 = Button()
        button2.Text = 'Exit'
        button2.AutoSize = True
        button2.Location = Point(375, lastpos+50)
        button2.Click += self.exit_button_clicked
        self.Controls.Add(button2)

    def on_enter(self, _sender, args):
        """
        Method invoked when a key is pressed within a textbox. It calls transform() if this key is enter
        :param: self the reference to the Form
        :param: _sender  ignore
        :param: args contains the pressed key event
        """
        if args.KeyCode == Keys.Enter:
            self.transform()

    def updatetbox_b(self, _sender, _event):
        """
        Method invoked when the beam angle slider is moved. Updates the text box and calls transform()
        :param: self the reference to the Form
        :param: _sender  ignore
        :param: _event ignore
        """
        self.tboxB.Text = str(self.tbB.Value)
        self.transform()

    def updatetbox_c(self, _sender, _event):
        """
        Method invoked when the couch angle slider is moved. Updates the text box and calls transform()
        :param: self the reference to the Form
        :param: _sender  ignore
        :param: _event ignore
        """
        self.tboxC.Text = str(self.tbC.Value)
        self.transform()

    def updatetbox_x(self, _sender, _event):
        """
        Method invoked when the x slider is moved. Updates the text box and calls transform()
        :param: self the reference to the Form
        :param: _sender  ignore
        :param: _event ignore
        """
        self.tboxX.Text = str(self.tbX.Value)
        self.transform()

    def updatetbox_y(self, _sender, _event):
        """
        Method invoked when the y slider is moved. Updates the text box and calls transform()
        :param: self the reference to the Form
        :param: _sender  ignore
        :param: _event ignore
        """
        self.tboxY.Text = str(self.tbY.Value)
        self.transform()

    def updatetbox_z(self, _sender, _event):
        """
        Method invoked when the z slider is moved. Updates the text box and calls transform()
        :param: self the reference to the Form
        :param: _sender  ignore
        :param: _event ignore
        """
        self.tboxZ.Text = str(self.tbZ.Value)
        self.transform()

    def updatetbox_e(self, _sender, _event):
        """
        Method invoked when the extraction slider is moved. Updates the text box and calls transform()
        :param: self the reference to the Form
        :param: _sender  ignore
        :param: _event ignore
        """
        self.tboxE.Text = str(self.tbE.Value)
        self.transform()

    def exit_button_clicked(self, _sender, _event):
        """
        Method invoked when the Exit button is clicked. It closes the form.
        :param: self the reference to the Form
        :param: _sender  ignore
        :param: _event ignore
        """
        self.Close()

    def flip_button_clicked(self, _sender, _event):
        """
        Method invoked when the Flip button is clicked. It toggles the flip boolean variable and calls the transform() function
        :param: self the reference to the Form
        :param: _sender  ignore
        :param: _event ignore
        """
        global flip
        flip = not flip
        self.transform()

    def apply_button_clicked(self, _sender, _event):
        """
        Method invoked when the Apply button is clicked. It calls the transform() function
        :param: self the reference to the Form
        :param: _sender  ignore
        :param: _event ignore
        """
        self.transform()

    def transform(self):
        """
        Slot function called whenever the Apply button is clicked, or when entered is clicked on text box,
        or when slider is moved so that text box is updated
        :param: self reference to the Form
        """
        # Get transformation from text box
        ba = self.tboxB.Text
        ca = self.tboxC.Text
        x = self.tboxX.Text
        y = self.tboxY.Text
        z = self.tboxZ.Text
        e = self.tboxE.Text if extraction else "0"

        # Sanity check that we are in the correct range
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
            z = str(int(self.tbZ.Minimum))
            self.tboxZ.Text = z
            ok = False
        if extraction:
            if e == "" or float(e) < self.tbE.Minimum:
                e = str(int(self.tbE.Minimum))
                self.tboxE.Text = e
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
        if extraction:
            if float(e) > self.tbE.Maximum:
                e = str(int(self.tbE.Maximum))
                self.tboxE.Text = e
                ok = False

        self.update_sliders()  # Update slider position

        # If input value was in correct interval, perform the transformation
        if ok:
            global gangle
            global cangle
            global bangle
            global tangle
            global oldgangle
            global oldcangle
            global oldbangle
            global oldtangle
            global cx
            global oldcx
            global cy
            global oldcy
            global cz
            global oldcz
            global se
            global oldse
            oldgangle = gangle
            oldcangle = cangle
            oldbangle = bangle
            oldtangle = tangle
            oldcx = cx
            oldcy = cy
            oldcz = cz
            oldse = se
            gangle = radians(float(ba))
            cangle = radians(float(ca))
            bangle = 0  # to be determined later
            tangle = 0  # to be determined later
            cx = float(x)
            cy = float(y)
            cz = float(z)
            se = float(e)
            # Convert couch deviation to cm (RayStation coordinates)
            cx /= 10.
            cy /= 10.
            cz /= 10.
            se /= 10.
            # Transform the models
            transform_models()

    def update_sliders(self):
        """
        Update the GUI sliders if after text box input finished.
        It has to be done without emitting new signal, to avoid an infinite loop
        :param: self reference to Form
        """
        # Get new values from text box
        newb = round(float(self.tboxB.Text))
        newc = round(float(self.tboxC.Text))
        newx = round(float(self.tboxX.Text))
        newy = round(float(self.tboxY.Text))
        newz = round(float(self.tboxZ.Text))
        newe = round(float(self.tboxE.Text)) if extraction else 0
        # If different from trackbar value, disconnect temporarily from slots and update the value
        if abs(newb-self.tbB.Value) > 0:
            self.tbB.ValueChanged -= self.updatetbox_b
            self.tbB.Value = newb
            self.tbB.ValueChanged += self.updatetbox_b
        if abs(newc - self.tbC.Value) > 0:
            self.tbC.ValueChanged -= self.updatetbox_c
            self.tbC.Value = newc
            self.tbC.ValueChanged += self.updatetbox_c
        if abs(newx - self.tbX.Value) > 0:
            self.tbX.ValueChanged -= self.updatetbox_x
            self.tbX.Value = newx
            self.tbX.ValueChanged += self.updatetbox_x
        if abs(newy - self.tbY.Value) > 0:
            self.tbY.ValueChanged -= self.updatetbox_y
            self.tbY.Value = newy
            self.tbY.ValueChanged += self.updatetbox_y
        if abs(newz - self.tbZ.Value) > 0:
            self.tbZ.ValueChanged -= self.updatetbox_z
            self.tbZ.Value = newz
            self.tbZ.ValueChanged += self.updatetbox_z
        if extraction:
            if abs(newe - self.tbE.Value) > 0:
                self.tbE.ValueChanged -= self.updatetbox_e
                self.tbE.Value = newe
                self.tbE.ValueChanged += self.updatetbox_e


def tune_models():
    """
    This function creates a GUI form with sliders for adjusting interactively the treatment head and couch position.
    Once the user presses exit, the form is closed and the imported 3D models are removed.
    """
    aform = TuneModelsForm()
    Application.Run(aform)
    # Form closed, remove now imported ROIs
    remove_models()


def transform_models():
    """
    This function transforms the imported 3D models to match a new gantry and couch angle, or couch position
    """
    # First, rotate the treatment head to the new angle
    moved = False
    if abs(cangle - oldcangle) > 0 or abs(gangle - oldgangle) > 0 or abs(se - oldse) > 0:
        for part in linac.parts:
            if part.active:
                roi_name = part.name
                b = -cs*(oldcangle+c0)
                b2 = cs*(cangle+c0)
                d = gs*(gangle - oldgangle)  # g0 cancels
                ey = gs*(se - oldse) if part.retractable else 0
                case.PatientModel.RegionsOfInterest[roi_name].TransformROI3D(Examination=examination, TransformationMatrix={
                    'M11': cos(d)*cos(b)*cos(b2)-sin(b)*sin(b2), 'M12': -sin(d)*cos(b2), 'M13': -cos(d)*sin(b)*cos(b2)-cos(b)*sin(b2), 'M14': iso.x-iso.x*(cos(d)*cos(b)*cos(b2)-sin(b)*sin(b2))+iso.y*sin(d)*cos(b2)+iso.z*(cos(d)*sin(b)*cos(b2)+cos(b)*sin(b2))    ,
                    'M21': sin(d)*cos(b)                       , 'M22':  cos(d)        , 'M23': -sin(d)*sin(b)                       , 'M24': iso.y-iso.x* sin(d)*cos(b)                        -iso.y*cos(d)        +iso.z* sin(d)*sin(b)                        - ey,
                    'M31': cos(d)*cos(b)*sin(b2)+sin(b)*cos(b2), 'M32': -sin(d)*sin(b2), 'M33': -cos(d)*sin(b)*sin(b2)+cos(b)*cos(b2), 'M34': iso.z-iso.x*(cos(d)*cos(b)*sin(b2)+sin(b)*cos(b2))+iso.y*sin(d)*sin(b2)+iso.z*(cos(d)*sin(b)*sin(b2)-cos(b)*cos(b2))    ,
                    'M41': 0                                   , 'M42':  0             , 'M43':  0                                   , 'M44': 1                                                                                                                      })
                moved = True
    # Then, move the couch to a new position
    if abs(cx - oldcx) > 0 or abs(cy - oldcy) or abs(cz-oldcz) > 0 or abs(cangle-oldcangle) > 0:
        for part in couch.parts:
            if part.active:
                roi_name = part.name
                dx = cx - oldcx
                dy = cy - oldcy
                dz = cz - oldcz
                if not part.moveX:
                    dx = 0
                if not part.moveY:
                    dy = 0
                if not part.moveZ:
                    dz = 0
                if not part.scissor:
                    if abs(dx) > 0 or abs(dy) > 0 or abs(dz) > 0:
                        case.PatientModel.RegionsOfInterest[roi_name].TransformROI3D(Examination=examination, TransformationMatrix={
                            'M11': 1, 'M12': 0, 'M13': 0, 'M14': dx,
                            'M21': 0, 'M22': 1, 'M23': 0, 'M24': dy,
                            'M31': 0, 'M32': 0, 'M33': 1, 'M34': dz,
                            'M41': 0, 'M42': 0, 'M43': 0, 'M44': 1})
                        moved = True

    if len(lsci) >= 2:  # scissor robot defined. Distances below are hard coded for the moment
        global bangle, tangle, oldbangle, oldtangle
        bs = 170  # cm
        lb = 120
        lt = 100
        rholim = lt + lb  # cm = 1.2 m plus 1 m
        bx = iso.x - bs*sin(cs*cangle)
        bz = iso.z + bs*cos(cs*cangle)
        oldbx = iso.x - bs*sin(cs*oldcangle)
        oldbz = iso.z + bs*cos(cs*oldcangle)
        tx = iso.x + dx0 + cx
        tz = iso.z + dz0 + cz
        xd = bx - tx
        zd = bz - tz
        rho = sqrt(xd*xd + zd*zd)
        failed = rho > rholim

        if failed:
            # no solution found
            # put the base opposite to ISO and the top towards it
            bangle = cangle + radians(180)
            tangle = cangle
        else:
            # solve SSS triangle https://www.mathsisfun.com/algebra/trig-solving-sss-triangles.html
            a = lt
            b = lb
            c = rho
            alpha = acos((b*b+c*c-a*a)/2/b/c)
            beta = acos((a*a+c*c-b*b)/2/c/a)
            delta = atan2(xd, zd)
            bangle = (delta + alpha)
            tangle = -(beta - delta)
            global flip
            if flip:
                bangle -= 2*alpha
                tangle += 2*beta

        for i, roi_name in enumerate(lsci):
            part = [p for p in couch.parts if p.name == roi_name][0]
            dx = cx - oldcx
            dy = cy - oldcy
            dz = cz - oldcz

            if i == 0:
                d = cs * (bangle - oldbangle)
            elif i == 1:
                d = cs * (tangle - oldtangle)
            else:
                d = cs * (cangle - oldcangle)

            if not part.moveX:
                dx = 0
            if not part.moveY:
                dy = 0
            if not part.moveZ:
                dz = 0

            if i == 0:
                rtpx = oldbx  # rotation point
                rtpz = oldbz  # rotation point
                dx = -bs*(sin(cs*cangle)-sin(cs*oldcangle))
                dz =  bs*(cos(cs*cangle)-cos(cs*oldcangle))
            elif i == 1:
                rtpx = iso.x + dx0 + oldcx
                rtpz = iso.z + dz0 + oldcz
            else:
                rtpx = iso.x
                rtpz = iso.z

            case.PatientModel.RegionsOfInterest[roi_name].TransformROI3D(Examination=examination, TransformationMatrix={
                'M11': cos(d), 'M12': 0, 'M13': -sin(d), 'M14': rtpx - rtpx*cos(d) + rtpz*sin(d) + dx,
                'M21': 0     , 'M22': 1, 'M23': 0      , 'M24': dy,
                'M31': sin(d), 'M32': 0, 'M33':  cos(d), 'M34': rtpz - rtpx*sin(d) - rtpz*cos(d) + dz,
                'M41': 0     , 'M42': 0, 'M43': 0      , 'M44': 1                                    })
            moved = True

    if moved:
        global colthread
        if colthread.IsAlive:
            colthread.Abort()
            colthread.Join()
        colthread = Thread(ThreadStart(detect_collision))
        colthread.Start()


def remove_models():
    """
    This function remove the ROIs created at the beginning of the script, to clean up everything upon script termination
    """
    for part in itertools.chain(linac.parts, couch.parts):
        if part.active:
            # delete ROI
            roi_name = part.name
            case.PatientModel.RegionsOfInterest[roi_name].DeleteRoi()


def detect_collision():
    """
    This function is used for automatic collision detection
    """
    try:
        print('Started')
    finally:
        Thread.Sleep(5000)
        print('Finished')


def main():
    """
    Start program
    """

    # Define your 3D models and machines available at your institution
    agility = Machine("Elekta Agility", "Y:\\ydrive\\STL parts\\Elekta Agility\\",
                      [Part("Nozzle", "RotatingHeads.stl", "Blue", True),
                       Part("Electron cone 15cm x 15cm", "ElectronApplicator.stl", "Blue", False),
                       Part("Flat panel left", "FlatPanelExtensionLeft.stl", "Blue", False),
                       Part("Flat panel bottom", "FlatPanelExtension.stl", "Blue", False),
                       Part("CBCT source", "CBCTsource.stl", "Blue", False)
                       ]
                      )
    proteus = Machine("iba Proteus", "Y:\\ydrive\\STL parts\\iba Proteus\\",
                      # [Part("Nozzle", "NozzleWithSmallSnout.stl", "Blue", True),
                      [Part("Nozzle", "Nozzle.stl", "Blue", True),
                       Part("Snout 18 cm", "Snout18cm.stl", "Blue", True, retractable=True),
                       Part("Range shifter 80mm 18cm", "RangeShifter80mm_18cm.stl", "Blue", False, retractable=True),
                       ]
                      )

    trilogy = Machine("Varian IX Trilogy", "Y:\\ydrive\\STL parts\\Varian IX Trilogy\\", [Part("Nozzle", "RotatingHeads.stl", "Blue", True)])
    truebeam = Machine("Varian TrueBeam", "Y:\\ydrive\\STL parts\\Varian TrueBeam\\", [Part("Nozzle", "RotatingHeads.stl", "Blue", True)])
    # For the couch model, you can specify which subparts are fixed in x, y, z coordinates, i.e. do not translate when moving the upper part of the couch
    evo = Machine("Hexapod Evo", "Y:\\ydrive\\STL parts\\Hexapod Evo\\",
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
    robot = Machine("Sciss Robot", "Y:\\ydrive\\STL parts\\Scissor Robot\\",
                    [Part("Couch", "Couch.stl", "Green", True),
                     Part("Scissor pedestal", "ScissorPedestal.stl", "Gray", True, False, True, False, True),
                     Part("Scissor base", "ScissorBase.stl", "Gray", True, False, True, False, True),
                     Part("Scissor top", "ScissorTop.stl", "White", True, True, True, True, True),
                     ]
                    )

    # Define the list of available treatment heads
    linacs = {agility.name: agility, proteus.name: proteus, trilogy.name: trilogy, truebeam.name: truebeam}
    # Define the list of available couches
    couches = {evo.name: evo, robot.name: robot}

    # GUI initialization. Some variables below are global

    # Select first which treatment head to use
    form = SelectListForm(linacs, "LINAC")
    Application.Run(form)
    global linac
    linac = linacs[form.name]

    # Select which subparts of the treatment head you want to draw
    pform = SelectPartsForm(linac)
    Application.Run(pform)

    # Select now which couch to use
    cform = SelectListForm(couches, "patient couch")
    Application.Run(cform)
    global couch
    couch = couches[cform.name]

    # Select which subparts of the couch you want to draw
    pcform = SelectPartsForm(couch)
    Application.Run(pcform)

    # Get now the currently opened RayStation patient
    global case
    case = get_current('Case')
    global examination
    examination = get_current('Examination')
    global structure_set
    structure_set = case.PatientModel.StructureSets[examination.Name]
    orientation = examination.PatientPosition

    # A rotation angle offset in degrees of the 3D model is needed to match the CT orientation depending on PatientPosition attribute
    # Also, a correction of the rotation direction is needed depending on the patient orientation.
    # TODO: change this and read instead the PatientOrientationMatrix from the DICOM CT
    gantry_angle_offset = {'HFS': 180, 'FFS': 180, 'HFP':   0, 'FFP':  0}
    couch_angle_offset =  {'HFS': 180, 'FFS':   0, 'HFP': 180, 'FFP':  0}
    gantry_direction =    {'HFS':  -1, 'FFS':  -1, 'HFP':  -1, 'FFP': -1}
    couch_direction =     {'HFS':  -1, 'FFS':  -1, 'HFP':   1, 'FFP':  1}
    # g0, c0 are the needed gantry angle and couch angle rotation of the 3D model to match this patient orientation
    # gs, cs are the rotation direction signs to be applied in order to match this patient orientation
    global g0, c0, gs, cs
    g0 = radians(gantry_angle_offset[orientation])
    c0 = radians(couch_angle_offset[orientation])
    gs = gantry_direction[orientation]
    cs = couch_direction[orientation]

    # Check if Isocenter has already been defined, if not, wait until defined, then continue
    poi_type = 'Isocenter'
    poi_lst = [r.Type for r in case.PatientModel.PointsOfInterest]
    while poi_type not in poi_lst:
        await_user_input('Please click OK and define an "'+poi_type+'" POI, then click on Play Script')
        poi_lst = [r.Type for r in case.PatientModel.PointsOfInterest]

    # If there are more than one isocenter, ask the user to confirm which one to use
    global iso
    if poi_lst.count(poi_type) > 1:
        isocenters = [r.Name for r in case.PatientModel.PointsOfInterest if r.Type == poi_type]
        print(isocenters)
        isolist = {isocenters[i]: i for i in range(0, len(isocenters))}
        print(isolist)
        isoform = SelectListForm(isolist, "Isocenter")
        Application.Run(isoform)
        iso = structure_set.PoiGeometries[isoform.name].Point
    else:
        iso = structure_set.PoiGeometries[poi_lst.index(poi_type)].Point

    # Create first model at angles g=0,c=0.
    # These below are global variables describing gantry angle (gangle), couch angle (cangle), couch position (cx,cy,cz), snout extration (se), and the old value before changing it.
    # tangle, bangle, lsci and flip are used just for scissor robot
    global gangle, cangle, bangle, tangle, oldgangle, oldcangle, oldbangle, oldtangle
    global cx, cy, cz, se, oldcx, oldcy, oldcz, oldse
    global flip
    global lsci
    gangle = 0
    cangle = 0
    bangle = 0
    tangle = 0
    oldgangle = 0
    oldcangle = 0
    oldbangle = 0
    oldtangle = 0
    cx = 0
    cy = 0
    cz = 0
    se = 0
    oldcx = 0
    oldcy = 0
    oldcz = 0
    oldse = 0
    flip = False
    lsci = []

    # Remove previous ROIs if already defined, e.g. if previous program instance crashed or script was stopped. This prevents an error later when importing.
    # User is asked for individual removal confirmation, just in case someone defined a clinical ROI with by chance the same name than your model.
    roi_lst = [r.Name for r in case.PatientModel.RegionsOfInterest]
    for part in itertools.chain(linac.parts, couch.parts):
        if part.active:
            roi_name = part.name
            if roi_name in roi_lst:
                await_user_input('Confirm deletion of preexisting ROI "' + roi_name + '" by clicking on Resume Script. Otherwise click Stop Script.')
                case.PatientModel.RegionsOfInterest[roi_name].DeleteRoi()

    # Create now treatment head ROIs and import STL models. Gantry and couch angle will be zero, and model will be centered at iso
    for part in linac.parts:
        if part.active:
            # create ROI
            roi_name = part.name
            roi_color = part.color
            roi_type = 'Support'
            file_name = linac.path + part.filename
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

    # Create now couch ROIs and import STL models. Couch will be centered at iso, but not moved.
    # Thus, it might be far away from the patient and has to be readjusted with the GUI sliders.
    for part in couch.parts:
        if part.active:
            # create ROI
            roi_name = part.name
            roi_color = part.color
            roi_type = 'Support'
            file_name = couch.path+part.filename
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

    # Check if a scissor robot is defined and store their part names in a list, being the first element the base, and the second element the top part,
    # and the third element the pedestal, if any
    auxlist = [p.name for p in couch.parts if p.scissor and p.active]
    lsci = []
    if len(auxlist) >= 2:
        lsci.append([pname for pname in auxlist if "base" in pname][0])
        lsci.append([pname for pname in auxlist if "top" in pname][0])
        lsci.append([pname for pname in auxlist if "pedestal" in pname][0])

    # Get list contoured couch ROIs here, ie. whose name contain couch (case insensitive)
    couch_lst = [r.Name for r in case.PatientModel.RegionsOfInterest if r.Type == 'Support' if re.search('couch', r.Name, re.IGNORECASE)]
    # Get list of couch STL 3D models, ie. whose name contain couch (case insensitive)
    couch_models = [c.name for c in couch.parts if c.active if re.search('couch', c.name, re.IGNORECASE)]

    # If there is a Couch ROI that someone contoured on the CT, recenter couch parts to match it approximately.
    # This is implemented by looking for the first occurrence ROI or model containing the substring couch.
    if len(couch_models) > 0:
        model_name = couch_models[0]
        geom = structure_set.RoiGeometries[model_name]
        if len(couch_lst) > 0:
            roi_name = couch_lst[0]
            geo = structure_set.RoiGeometries[roi_name]

            # Get center of the contoured couch
            rb = geo.GetBoundingBox()
            rx = (rb[1].x + rb[0].x) / 2
            ry = (rb[1].y + rb[0].y) / 2
            rz = rb[0].z
            # Get center of the 3D model couch
            mb = geom.GetBoundingBox()
            mx = (mb[1].x + mb[0].x) / 2
            my = (mb[1].y + mb[0].y) / 2
            mz = mb[0].z

            # Move all couch parts to close the difference
            global dx0
            global dz0
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
                    if not part.scissor:
                        case.PatientModel.RegionsOfInterest[roi_name].TransformROI3D(Examination=examination, TransformationMatrix={
                            'M11': 1, 'M12': 0, 'M13': 0, 'M14': dx,
                            'M21': 0, 'M22': 1, 'M23': 0, 'M24': dy,
                            'M31': 0, 'M32': 0, 'M33': 1, 'M34': dz,
                            'M41': 0, 'M42': 0, 'M43': 0, 'M44': 1})
                    else:
                        case.PatientModel.RegionsOfInterest[roi_name].TransformROI3D(Examination=examination, TransformationMatrix={
                            'M11': 1, 'M12': 0, 'M13': 0, 'M14': dx,
                            'M21': 0, 'M22': 1, 'M23': 0, 'M24': dy,
                            'M31': 0, 'M32': 0, 'M33': 1, 'M34': dz,
                            'M41': 0, 'M42': 0, 'M43': 0, 'M44': 1})

    global extraction
    extraction = any([part.retractable for part in linac.parts if part.active])

    # Collision detection thread
    global colthread
    colthread = Thread(ThreadStart(detect_collision))

    # Tuning form thread
    thread = Thread(ThreadStart(tune_models))
    thread.Start()
    thread.Join()
    if colthread.IsAlive:
        colthread.Abort()
        colthread.Join()

if __name__ == '__main__':
    main()
