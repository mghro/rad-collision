# RadCollision

RadCollision is an open-source platform for designed to aid treatment planners in the choice of beam incidence angles that do not collide with the patient or couch. Currently, only the RayStation TPS is supported.
It loads a 3D model of a LINAC into RayStation, that can be rotated interactively to assess colission risk.
It is written in the native IronPython scripting language of the RayStation interface.

Licensing
---------

This software platform is open-source. You can use, modify and contribute to it. Original attribution/citation shall be granted to:

Hueso-González et al. - Biomed Phys Eng Express 2020, "An open-source platform for interactive collision prevention in photon and particle beam therapy treatment planning".


Requirements
------------

- Raystation 8B or higher
- 3D model of your nozzle and couch (optional) as STL files


How to use
----------

- Create a new script in RayStation and paste the contents of collision_detection.py
- Modify the path where the STL files are stored, to point at your own drive


3D model format
---------------

- The file type should be STL
- The model origin shall be in the isocenter
- The perspective should match that of an observer standing in front of the gantry
- More information is given in the header of collision_detection.py

Known issues
------------

- It does not handle yet custom PatientOrientations
- RayStation (9A) crashes if you import STL geometries that are not closed surfaces
- RayStation (9A) crashes when calculating DSC for very complex geometries with holes


Created and documented by
-------------------------

Fernando Hueso-Gonzalez - fhuesogonzalez@mgh.harvard.edu
Massachusetts General Hospital and Harvard Medical School

Acknowledgment
--------------

Erik Vidholm (RaySearch)
