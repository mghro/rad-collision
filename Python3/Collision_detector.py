# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 19:15:29 2021

@author: kpiqu
"""

import numpy as np
import math
import stl
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot
from math import cos, sin, radians
import tkinter as tk
from tkinter import *

import vtk
import argparse
import itertools
import vedo


#The values of the starting set-up:
anggantry = radians(0)
angcouch = radians(0)
posisocent =[0,0,0]
posisocent = np.array(posisocent)

#Here we define the functions that are going to be used.

def translate(self, translation):
    '''
    Translate the mesh in the three directions
    :param numpy.array translation: Translation vector (x, y, z)
    '''
    assert len(translation) == 3, "Translation vector must be of length 3"
    self.x += translation[0]
    self.y += translation[1]
    self.z += translation[2]

def Matriz_Rotacion (axis, angle):
    """
    Parameters
    ----------
    axis : string
        String varible that especifies the axis in which the rotation matrix is going to be produced.
    angle : float
        The angle that is going to be rotated.

    Raises
    ------
    RuntimeError
        DESCRIPTION.

    Returns
    -------
    matrix : np.array like matriz
        DESCRIPTION.

    """
    
    matrix = []
    if 'x' == axis:
        matrix = [[1,0,0],[0,cos(angle), sin(angle)],[0, -sin(angle),cos(angle)]]
    elif 'y' == axis:
        matrix = [[cos(angle),0, -sin(angle)],[0,1,0],[sin(angle),0 ,cos(angle)]]
    elif 'z' == axis:
        matrix = [[cos(angle), sin(angle),0],[-sin(angle),cos(angle), 0],[0,0,1]]
    else:
        raise RuntimeError('Unknown axis %r, expected x, y or z' % axis)
    
    matrix = np.array(matrix)
    return matrix

def rotate_using_matrix(self, rotation_matrix, point=None):
    '''
    Rotate using a given rotation matrix and optional rotation point

    Note that this rotation produces clockwise rotations for positive
    angles which is arguably incorrect but will remain for legacy reasons.
    For more details, read here:
    https://github.com/WoLpH/numpy-stl/issues/166
    '''
    identity = np.identity(rotation_matrix.shape[0])
    # No need to rotate if there is no actual rotation
    if not rotation_matrix.any() or (identity == rotation_matrix).all():
        return

    if isinstance(point, (numpy.ndarray, list, tuple)) and len(point) == 3:
        point = numpy.asarray(point)
    elif point is None:
        point = numpy.array([0, 0, 0])
    elif isinstance(point, (int, float)):
        point = numpy.asarray([point] * 3)
    else:
        raise TypeError('Incorrect type for point', point)

    def _rotate(matrix):
        if point.any():
        # Translate while rotating
            return (matrix - point).dot(rotation_matrix) + point
        else:
        # Simply apply the rotation
            return matrix.dot(rotation_matrix)

        # Rotate the normals
        self.normals[:] = _rotate(self.normals[:])

        # Rotate the vectors
        for i in range(3):
            self.vectors[:, i] = _rotate(self.vectors[:, i])

def find_mins_maxs(obj):
    """
    Parameters
    ----------
    obj : mesh
        Mesh file that is going to be analized.

    Returns
    -------
    minx : float
        DESCRIPTION.
    maxx : float
        DESCRIPTION.
    miny : float
        DESCRIPTION.
    maxy : float
        DESCRIPTION.
    minz : float
        DESCRIPTION.
    maxz : float
        DESCRIPTION.

    """
    minx = obj.x.min()
    maxx = obj.x.max()
    miny = obj.y.min()
    maxy = obj.y.max()
    minz = obj.z.min()
    maxz = obj.z.max()
    return minx, maxx, miny, maxy, minz, maxz

def mesh_to_vtkPolydata(obj):
    """
    Parameters
    ----------
    obj : mesh
        Mesh type variable with the information of the solid.
    
    Returns
    -------
    vpoly : VtkPolyDataObject
        VtkPolyDataObject type variable with the information of the solid..

    """
    verts = list(itertools.chain(*(obj.vectors)))
    faces = [[i*3, i*3+1, i*3+2] for i in range(len(verts)//3)]
    vpoly = vedo.Mesh([verts, faces]).clean().polydata()
    return vpoly

def Choque(a,b):
    """
    Parameters
    ----------
    a : mesh
        Mesh type variable with the information of the solid which interaction will be calculated.
    
    b : mesh
        Mesh type variable with the information of the solid which interaction will be calculated..
       
    Returns
    -------
    Boolean object that especifies if there is a collision between those two stl files.

    """
    a = mesh_to_vtkPolydata(a)
    b = mesh_to_vtkPolydata(b)
                
    #Calculamos si se realiza la colisión
    collide = vtk.vtkCollisionDetectionFilter()
    collide.SetInputData(0, a)
    collide.SetInputData(1, b)
    collide.SetTransform(0, vtk.vtkTransform())
    collide.SetMatrix(1, vtk.vtkMatrix4x4())
    collide.SetBoxTolerance(0.0)
    collide.SetCellTolerance(0.0)
    collide.SetNumberOfCellsPerNode(2)
    collide.SetCollisionModeToFirstContact()
    collide.GenerateScalarsOn()
    collide.Update()
      
    if collide.GetNumberOfContacts() > 0:
        return True
    else:
        return False

def get_program_parameters():
    description = 'Collision detection.'
    epilogue = '''

    '''
    parser = argparse.ArgumentParser(description='Process some integers.')

    parser.add_argument('--fileGantry', dest='fileGantry', type=str, required=True)
    parser.add_argument('--fileCouch', dest='fileCouch', type=str, required=True)
    parser.add_argument('--fileBody', dest='fileBody', type=str, required=True)
    parser.add_argument('--RotPat', dest='RotPat',default=False, type=bool)
    args = parser.parse_args()


    return args

# =============================================================================
# def Distancia(a,b):
    
    
    
   # a = mesh_to_vtkPolydata(a)
   # b = mesh_to_vtkPolydata(b)
    
   # distanceFilter = vtk.vtkDistancePolyDataFilter()
   # distanceFilter.SetInputData(0,a)
   # distanceFilter.SetInputData(1,b)
   # distanceFilter.Update()
    
   # rango= distanceFilter.GetOutput().GetPointData().GetScalars().GetRange() #Se supone que esta función calcula la mínima y la máxima distancia, he comprobado y creo que está bien
   # mindist = rango[0]

   # return mindist
# =============================================================================



def plot_something(gantry,couch,body,wix,wiy,wiz,wpx,wpy,w1,w2,figure):
    """
    

    Parameters
    ----------
    gantry :  Mesh
        DESCRIPTION.
    couch : Mesh
        DESCRIPTION.
    body : Mesh
        DESCRIPTION.
    wix : Scale
        Scale object of tkinter module.
    wiy : Scale
        Scale object of tkinter module.
    wiz : Scale
        Scale object of tkinter module.
    wpx : Scale
        Scale object of tkinter module.
    wpy : Scale
        Scale object of tkinter module.    
    w1 : Scale
        Scale object of tkinter module.
    w2 : Scale
        Scale object of tkinter module.
    figure : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

        
    #First, we locate the isocenter.Then we move the patient and the couch to the correct position. 
    veciso = [wix.get(),wiy.get(),wiz.get()]
    veciso = np.array(veciso)
    vecpac = [wpx.get(),wpy.get(),0]
    vecpac = np.array(vecpac)
    
    couch.translate(-veciso -vecpac)
    body.translate(-veciso)

    cvolume, ccog, cinertia = couch.get_mass_properties()
    
    #Now, we allow the couch-body sistem to move
    
    CMatriz = Matriz_Rotacion('z',math.radians(w2.get()))
    CMatriz = np.array(CMatriz)
    
    global poscouch
    
    couch.rotate_using_matrix(CMatriz,ccog)
    body.rotate_using_matrix(CMatriz,ccog -poscouch)    
    CMatrizInv= np.linalg.inv(CMatriz)
    
    #As the isocenter moves the gantry has to translate
    newveciso = np.dot(CMatriz,veciso)
    gantry.translate(newveciso-veciso) 
    
    #Once the gantry has been translated above the isocenter, now we allow it to move.
    GMatriz = Matriz_Rotacion('y',math.radians(w1.get()))
    GMatriz = np.array(GMatriz)
    gantry.rotate_using_matrix(GMatriz)
    GMatrizInv= np.linalg.inv(GMatriz) 
    
    
    #We evaluate if there is a collision
    HayChoque1 = Choque(gantry,couch)
    HayChoque2 = Choque(gantry,body)
    HayChoque3 = Choque(body,couch)
    HayChoque = HayChoque1 or HayChoque2 or HayChoque3
    
    if HayChoque==False:
        print("No collision in this configuration")
        
        global axes
        if axes.collections:
            axes.collections[-1].remove()
            axes.collections[-1].remove()
            axes.collections[-1].remove()  
        
        
        axes.add_collection3d(mplot3d.art3d.Poly3DCollection(gantry.vectors))
        axes.add_collection3d(mplot3d.art3d.Poly3DCollection(couch.vectors,facecolors='green'))
        axes.add_collection3d(mplot3d.art3d.Poly3DCollection(body.vectors, facecolors='red'))
    
        figure.canvas.draw()
        
        #Calculamos la distancia mínima
        #Distancia1 = Distancia(gantry,couch)
        #Distancia2 = Distancia(gantry,body)
        #Distancia3 = Distancia(body,couch)
        #print(f'La distancia mínima entre gantry y couch es: {Distancia1:,.0f} mm')
        #print(f'La distancia mínima entre gantry y paciente es: {Distancia2:,.0f} mm')
        
        
    else:
        print("This set of conditions results in a collision. Insert other values")
        
    
    
        
    
    #Once we have finished evaluating the configuration, we reset the sistem to his original configuration
    #befour the user press again the 'Show' Button. 

    gantry.rotate_using_matrix(GMatrizInv)
    gantry.translate(-newveciso+veciso)
    couch.rotate_using_matrix(CMatrizInv, ccog)    
    body.rotate_using_matrix(CMatrizInv, ccog -poscouch) 
    couch.translate(+veciso)
    body.translate(+veciso)

def main():   
    
    #Lets read the parameters from the command line
    args = get_program_parameters()

    # Load the STL files and add the vectors to the plot
    gantry = mesh.Mesh.from_file(args.fileGantry)
    couch = mesh.Mesh.from_file(args.fileCouch)
    body = mesh.Mesh.from_file(args.fileBody)
    

    # Move STL file from RayStation coordinate system to IEC one.
    MatrizPrueba = Matriz_Rotacion('x',radians(90))
    gantry.rotate_using_matrix(MatrizPrueba)
    couch.rotate_using_matrix(MatrizPrueba)
    body.rotate_using_matrix(MatrizPrueba)
    
    if args.RotPat == True:
        MatrizPrueba = Matriz_Rotacion('y',radians(180))
        body.rotate_using_matrix(MatrizPrueba)
    
    #We evaluate the body and we locate the couch underneath it.
    bvolume, bcog, binertia = body.get_mass_properties()
    body.translate(-posisocent)
    
    bminx, bmaxx, bminy, bmaxy, bminz, bmaxz = find_mins_maxs(body)
    bw1 = bmaxx - bminx
    bl1 = bmaxy - bminy
    bh1 = bmaxz - bminz
    
    global poscouch
    poscouch = [0,0,-bh1/1.5] 
    
    cvolume, ccog, cinertia = couch.get_mass_properties()
    couch.translate(-ccog + poscouch -posisocent)
    cvolume, ccog, cinertia = couch.get_mass_properties()
    couch.rotate_using_matrix(Matriz_Rotacion('z',angcouch),ccog)
    
    MatrizPrueba = Matriz_Rotacion('y', anggantry)
    gantry.rotate_using_matrix(MatrizPrueba,[0,0,0])
    
    # Create a new plot
    figure = pyplot.figure()
    global axes
    axes = mplot3d.Axes3D(figure)
    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(gantry.vectors))
    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(couch.vectors,facecolors='green'))
    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(body.vectors, facecolors='red'))
    
    # Auto scale to the mesh size
    scale1 = body.points.flatten()
    scale2 = gantry.points.flatten()
    scale3 = couch.points.flatten()
    scale = np.concatenate([scale1,scale2,scale3])
    axes.auto_scale_xyz(scale, scale, scale)
    axes.set_xlabel('$X$')
    axes.set_ylabel('$Y$')
    axes.set_zlabel('$Z$')        
    
    # Show the plot to the screen
    pyplot.show()
    
    #Activate the Slider interface
    master = tk.Tk()
    master.geometry("600x800")
    
    my_label1 = tk.Label(master,text="Gantry angle").place(x=250, y = 0)
    w1 = tk.Scale(master, from_=0, to=360, tickinterval=90)
    w1.place(x = 250, y = 30)
    
    my_label2 = tk.Label(master,text="Couch angle").place(x = 250, y = 150)
    w2 = tk.Scale(master, from_=0, to=180, tickinterval=45)
    w2.place(x = 250, y = 180)
    
    my_label3 = tk.Label(master,text="Isocenter position").place(x = 250, y = 300)
    
    my_labelix = tk.Label(master,text="X").place(x = 50, y = 330)
    wix = tk.Scale(master,  from_=-500, to=500, tickinterval=250)
    wix.place(x = 50, y = 360)
    
    my_labeliy = tk.Label(master,text="Y").place(x = 250, y = 330)
    wiy = tk.Scale(master, from_=-1500, to=1500, tickinterval=750)
    wiy.place(x = 250, y = 360)
    
    my_labeliz = tk.Label(master,text="Z").place(x = 450, y = 330)
    wiz = tk.Scale(master, from_=-500, to=500, tickinterval=250)
    wiz.place(x = 450, y = 360)
        
    my_label4 = tk.Label(master,text="Patient position").place(x = 250, y = 520)
    
    my_labelpx = tk.Label(master,text="X").place(x = 50, y = 550)
    wpx = tk.Scale(master, from_=-150, to=150, tickinterval=75)
    wpx.place(x = 50, y = 580)
    
    my_labelpy = tk.Label(master,text="Y").place(x = 250, y = 550)
    wpy = tk.Scale(master, from_=-150, to=150, tickinterval=75)
    wpy.place(x = 250, y = 580)    
    
    boton1= tk.Button(master, text='Show', command= lambda: plot_something(gantry,couch,body,wix,wiy,wpx,wpy,wiz,w1,w2,figure))
    boton1.place(x = 250, y = 720)
    
    tk.mainloop()
    


if __name__ == "__main__":
    main()
