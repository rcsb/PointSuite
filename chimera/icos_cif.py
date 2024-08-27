# -----------------------------------------------------------------------------
# Covering or iconifying the Chimera window while it is rendering images
# will corrupt the image.  
#
# based on code written for Chimera version 1.2199 by Tom Goddard, 3/27/06,
# (goddard@cgl.ucsf.edu). -- C. Lawson August 2006
# 
# updated 15 June 2011 with changes from Batsal Devkota
# updated in 2012 by Ezra Peisach to handle split entries
# in development: icos virus capsid from CIF file in place of PDB
# this almost works because chimera reads in the integer id matrices
# from _pdbx_struct_oper_list and assumes they correspond to BIOMT.
# in ~90 percent of cases this is true!
#

from os      import environ
from os.path import join, basename
from chimera import viewer, runCommand, openModels
import chimera

cif_files = environ["CHIMERAPDB"]

red = (1,0,0,1)         # Red, green, blue, opacity
green = (0,1,0,1)
blue = (0,0,1,1)
cyan = (0,1,1,1)
yellow=(1,1,0,1)
magenta=(1,0,1,1)
grey=(1,1,1,1)
resolution = 5          # Surface resolution, angstroms.
dens = 0.02             # Controls surface size.  Atoms / A^3.
dens_ca_only = 0.002    # Controls surface size for CA chains.
supersample = 3         # Capture images N times bigger for better quality.
padding = .05           # Fraction of window size for padding at edge.

import MultiScale
d = MultiScale.multiscale_model_dialog(create = True)
d.surface_resolution.set('%.4g' % resolution)
d.density_threshold.set('%.6g' % dens)
d.density_threshold_ca_only.set('%.6g' % dens_ca_only)
d.multimer_type.set(d.multimer_biounit)
viewer.windowSize = (512, 512)  # Set Chimera window size
runCommand('set bg_color white')


def make_images(cif_files):
     # open structures
     molecules = []
     for i in cif_files:
          molecule = openModels.open(i, 'CIF/mmCIF')[0]
          molecules.append(molecule)

     #delete solvent
     runCommand('delete solvent')

     # Create multiscale surfaces
     d.make_multimers(molecules)

     #--------begin new commands added 2011----------------- 
     # Show only surfaces, even for original copy of molecule.
     d.select_all_chains_cb()
     d.show_only_surface_cb()
     d.clear_selection_cb()

     #Turn on silhouette edges 
     runCommand('set silhouette')

     #Use high quality lighting, requires Chimera v 1.4 or higher
     runCommand('set light_quality glossy')
     #--------end new commands added 2011----------------- 

     # turn model (better view for helical viruses)
     #runCommand('turn x -90')

     # Set camera so virus fills window.
     set_camera_parameters(padding)

     # Capture and save image
     from chimera.printer import saveImage
     for i in cif_files:
          image_path = join(i + '.jpg')
          saveImage(image_path, format='jpeg', supersample = supersample)

# ---------------------------------------------------------------------------
# This adjusts camera so virus capsid fills window.
#
def set_camera_parameters(padding):

    have_bbox, bbox = chimera.openModels.bbox()
    if not have_bbox:
        return False

    z1 = bbox.llf.z
    z2 = bbox.urb.z             # z2 > z1
    zsize = z2 - z1
    xsize = bbox.urb.x - bbox.llf.x
    ysize = bbox.urb.y - bbox.llf.y

    v = chimera.viewer
    c = v.camera
# new
    s = .5 * max(xsize, ysize)
    s += padding * s
    v.viewSize = s
#end new
    v.viewSize = .5 * max(xsize, ysize)
    v.scaleFactor = 1
    c.center = bbox.center().data()
    c.viewDistance = zsize
    c.nearFar = (z2 + padding*zsize, z1 - padding*zsize)

    return True
# -----------------------------------------------------------------------------
# Color certain matrix subunits and chains of a multiscale model.
#

def color_multiscale(rgba, model_num, matrix_nums, chain_ids = None):


  import MultiScale
  d = MultiScale.multiscale_model_dialog(True)
  for cp in d.chain_pieces():
    g = cp.surface_group
    if g.model().id == model_num:
      mnum, chain_id = g.oslName.split('.')
      if chain_ids == None or chain_id in chain_ids:
        if int(mnum) in matrix_nums:
          cp.set_color(rgba)


# ----------------------------------------------------------------------------
#
# This function call makes the image
#
make_images(cif_files.split())

#Exit Chimera
runCommand('stop noask')

