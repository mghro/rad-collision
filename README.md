[//]: # (RadCollision)

<img src="logo.svg" alt="Logo" width="200"/>

RadCollision is an open-source platform designed to aid treatment planners in the choice of beam incidence angles that do not collide with the patient or couch.
It loads a 3D model of a LINAC and or couch into the TPS, that can be rotated interactively to assess collision risk. One can also load in the TPS a 3D scan of the patient done e.g. with a smartphone, if full patient geometry is needed.

Here are some screenshots of the capabilities and a news release:
* https://www.youtube.com/watch?v=TlFGaw2NK6Q
* https://www.youtube.com/watch?v=6osJI_xftAE
* https://physicsworld.com/a/open-source-software-detects-potential-collisions-in-radiotherapy-plans/

The (only available) interface layer is written in the native IronPython scripting language of the RayStation interface. Extension to other TPS is envisioned. 

Alternatives
------------

Alternative solutions are described here:
* References of https://doi.org/10.1088/2057-1976/aba442
* SlicerRT: http://perk.cs.queensu.ca/sites/perkd7.cs.queensu.ca/files/Suriyakumar2017a.pdf
* Eclipse: https://doi.org/10.1002/acm2.12673
* Pinnacle: https://doi.org/10.1002/acm2.12915 https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.4924827
* Visual Basic: https://doi.org/10.1016/j.canrad.2020.01.003
* Visual C#.Net: https://doi.org/10.1002/acm2.12963
* MATLAB: https://doi.org/10.1002/acm2.12998
* Unity 3D: https://doi.org/10.3389/fonc.2021.617007
* Varian: https://doi.org/10.1002/acm2.13496

Licensing
---------

This software platform is open-source, and is designed to work / be embedded with commercial TPS through their native scripting interface. You can use, modify and contribute to it. Original attribution/citation shall be granted to:

F Hueso-González et al 2020 - Biomed. Phys. Eng. Express 6 055013, "An open-source platform for interactive collision prevention in photon and particle beam therapy treatment planning". https://doi.org/10.1088/2057-1976/aba442 https://arxiv.org/abs/2007.05248

For a fully open-source solution for collision detection including the TPS planning, please refer to https://github.com/SlicerRt/SlicerRT
http://perk.cs.queensu.ca/sites/perkd7.cs.queensu.ca/files/Suriyakumar2017a.pdf


Requirements
------------

- Raystation 8B or higher
- 3D model of your nozzle and couch (optional) as STL files

For first attempts, you can use the open-source files stored in SlicerRT:
https://github.com/SlicerRt/SlicerRT/tree/master/RoomsEyeView/TreatmentMachineModels


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

- See https://github.com/mghro/rad-collision/issues
- It does not handle yet custom PatientOrientations
- RayStation (10A) crashes if you import STL geometries that are not closed surfaces
- RayStation (10A) crashes when calculating DSC for very complex geometries with holes


Created and documented by
-------------------------

Fernando Hueso-González - fhuesogonzalez@mgh.harvard.edu
Massachusetts General Hospital and Harvard Medical School

Acknowledgment
--------------

Thanks to everyone reporting issues. See also Acknowledgment section on https://doi.org/10.1088/2057-1976/aba442 https://arxiv.org/abs/2007.05248
