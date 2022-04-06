# -*- coding: utf-8 -*-
"""
Created on Sat Aug 21 03:23:35 2021

@author: kpiqu
"""
import numpy as np
import math
import stl
from stl import mesh
from matplotlib import pyplot
from math import cos, sin, radians


import vtk
import argparse #Faltaría ponerlo que tb se pueda pasar argparse por aquí
import itertools
import vedo


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

def get_program_parameters():
    import argparse
    description = 'Collision detection.'
    epilogue = '''

    '''
    parser = argparse.ArgumentParser(description='Process some integers.')

    parser.add_argument('--fileGantry', dest='fileGantry', type=str)
    parser.add_argument('--fileCouch', dest='fileCouch', type=str)
    parser.add_argument('--fileBody', dest='fileBody', type=str)
    parser.add_argument('--RotPat', dest='RotPat',default=False, type=bool)
    args = parser.parse_args()

    return args

def Choque(a,b):
    """
    

    Parameters
    ----------
    a : mesh
        Mesh variable that is going to be analized.
    
    b : mesh
        Mesh variable that is going to be analized.
       
    Returns
    -------
    Boolean object that especifies if there is a collision between those two stl files.

    """
                
    vpoly = mesh_to_vtkPolydata(a)
    vpoly2 = mesh_to_vtkPolydata(b)
    
    collide = vtk.vtkCollisionDetectionFilter()
    collide.SetInputData(0, vpoly)
    collide.SetInputData(1, vpoly2)
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




def main():
       
    #Lets read the parameters from the command line
    args = get_program_parameters()

    # Load the STL files and add the vectors to the plot
    gantry = mesh.Mesh.from_file(args.fileGantry)
    couch = mesh.Mesh.from_file(args.fileCouch)
    body = mesh.Mesh.from_file(args.fileBody)
    
        
    #The values of the starting set-up:
    anggantry = radians(0)
    angcouch = radians(0)
    veciso =[0,800,0]
    veciso = np.array(veciso)
    
        
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
    body.translate(-veciso)
    
    bminx, bmaxx, bminy, bmaxy, bminz, bmaxz = find_mins_maxs(body)
    bw1 = bmaxx - bminx
    bl1 = bmaxy - bminy
    bh1 = bmaxz - bminz
    
    global poscouch
    poscouch = [0,0,-bh1/1.5] 
    
    cvolume, ccog, cinertia = couch.get_mass_properties()
    couch.translate(-ccog + poscouch-veciso)
    cvolume, ccog, cinertia = couch.get_mass_properties()#we locate the final center of gravity of the couch
    couch.rotate_using_matrix(Matriz_Rotacion('z',angcouch),ccog) 
    
    
    RegVali=[]
    RegNoVali = []
    cont = 0
    for i in range(36):#This will count the gantry angle (from 0 to 359)
        for j in range(36):#This will count the couch angle (from 0 to 359) 
            
            PAngG=i*10
            print(PAngG)
            PAngC= j*10
            
            #We allow the couch-body to rotate
    
            CMatriz = Matriz_Rotacion('z',math.radians(PAngC))
            CMatriz = np.array(CMatriz)
            
            couch.rotate_using_matrix(CMatriz,ccog)
            body.rotate_using_matrix(CMatriz,ccog -poscouch)    
            CMatrizInv= np.linalg.inv(CMatriz) 
            
                    
            #As the isocenter moves the gantry has to translate
            newveciso = np.dot(CMatriz,veciso)
            gantry.translate(newveciso-veciso) 
            #Once the gantry has been translated above the isocenter, now we allow it to move.
            GMatriz = Matriz_Rotacion('y',math.radians(PAngG))
            GMatriz = np.array(GMatriz)
            gantry.rotate_using_matrix(GMatriz)
            GMatrizInv= np.linalg.inv(GMatriz) 
            
            #We evaluate if there is a collision
            HayChoque1 = Choque(gantry,couch)
            HayChoque2 = Choque(gantry,body)
            HayChoque3 = Choque(body,couch)
            HayChoque = HayChoque1 or HayChoque2 or HayChoque3
            print(HayChoque)
            
    
            if HayChoque== False:
                RegVali.append([PAngG,PAngC])
            else:
                RegNoVali.append([PAngG,PAngC])
                
            cont = cont + 1
            
            #Once we have finished evaluating the configuration, we reset the sistem to his original configuration
            
            gantry.rotate_using_matrix(GMatrizInv)
            gantry.translate(-newveciso+veciso)
            couch.rotate_using_matrix(CMatrizInv, ccog)    
            body.rotate_using_matrix(CMatrizInv, ccog -poscouch) 
    
    
    #Now we calculate the percentajes of configurations
    PorcentajeValido = len(RegVali)/cont *100
    PorcentajeNoValido = len(RegNoVali)/cont*100
    
    
    #Create a new plot
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    
    data = np.array(RegVali)
    xval, yval = data.T
    datano = np.array(RegNoVali)
    xnoval, ynoval = datano.T
    
    pyplot.scatter(xval,yval,c='b', marker = 's',linewidth=0.01 , label=f'No hay colisión {PorcentajeValido:,.0f} %')
    pyplot.scatter(xnoval,ynoval,c='r', marker = 's',linewidth=0.01 , label=f'Hay colisión {PorcentajeNoValido:,.0f} %')
    
    ratio = 1.0
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)
    
    
    pyplot.xlabel("Ángulo del gantry")
    pyplot.ylabel("Ángulo de la camilla")
    pyplot.legend(bbox_to_anchor=(0.4, 1.15));
    pyplot.show()


if __name__ == "__main__":
    main()