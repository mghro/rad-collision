[//]: # (RadCollision-RayStation)

This folder contains one flavor (or interface layer) of RadCollision. It supports both IronPython 2.7 as well as CPython 2.7 and 3.6, the latter being the only option starting from RayStation version 12B. For other languages, consult https://github.com/mghro/rad-collision

Licensing
---------

Please refer to https://github.com/mghro/rad-collision/blob/main/README.md and https://github.com/mghro/rad-collision/blob/main/LICENSE

Citation
--------

F Hueso-González et al 2020 - Biomed. Phys. Eng. Express 6 055013, "An open-source platform for interactive collision prevention in photon and particle beam therapy treatment planning". https://doi.org/10.1088/2057-1976/aba442 https://arxiv.org/abs/2007.05248

Requirements
------------

- Raystation version 8B or higher
- 3D model of your nozzle and couch (optional) as STL files

For first attempts, you can use the open-source STL files stored in this [PR](https://github.com/mghro/rad-collision/issues/21#issuecomment-1073840985) or in [https://github.com/SlicerRt/SlicerRT/tree/master/RoomsEyeView/TreatmentMachineModels](SlicerRT).

How to use
----------

- Create a new script in RayStation and paste the contents of collision_detection.py
- Modify the path where the STL files are stored, to point at your own drive


[![screenshot](screenshot.png)](screenshot.png)


3D model format
---------------

- The file type should be STL
- The model origin shall be in the isocenter
- The perspective should match that of an observer standing in front of the gantry
- More information is given in the header of collision_detection.py

Known issues
------------

- See https://github.com/mghro/rad-collision/issues
- It does not handle yet custom PatientOrientations
- RayStation (10A) crashes if you import STL geometries that are not closed surfaces
- RayStation (10A) crashes when calculating DSC for very complex geometries with holes

Authors
-------

- F. Hueso-González, MGH, now at IFIC (CSIC-UV)
- J. Morin, CISSS-CA
