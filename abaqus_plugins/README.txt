================================================================================
								Abaqus plugin
================================================================================

1. The folder 'WovenComposite' containes the abaqus plugin for woven composite 
   material assignment.
2. Copy the folder to the abaqus_plugin directory/current work directory.
3. The plugin requires the 'Eshebly.so' python module to work.
4. If any error realted to import of Eshelby module is encountered, try again by 
   compiling the Eshebly module using the make file by specifing the current 
   version of Abaqus being used.
5. If any error related to 'rsgTmpDir' is encountered, go to abaqus_plugin folder
   and remove the folder 'rsgTmpDir' and try again. 
