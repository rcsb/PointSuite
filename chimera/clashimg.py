# -----------------------------------------------------------------------------
# Script to make a crystal contact image for a PDB file.
# C Lawson, Sept 1 2006
#
# modification of a script written by Tom Goddard on Aug 8, 2006.
#

from os      import environ
from os.path import join, basename
from chimera import viewer, runCommand, openModels

pdb_path = environ["CHIMERAPDB"]
                                                                                

viewer.windowSize = (512, 512)  # Set Chimera window size
runCommand('set bg_color white')

image_path = join(pdb_path + '.png')
molecule = openModels.open (pdb_path, 'PDB')[0]
runCommand('crystalcontacts #0 1.5')
#runCommand('copy file %s png' % image_path)
supersample = 3
from chimera.printer import saveImage
saveImage(image_path, format = 'png', supersample = supersample)
