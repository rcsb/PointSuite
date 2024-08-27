#!/bin/csh -f
#
# this version edited by Ezra to handle split entries
##########################################################################
# home of chimera scripts relative to demo:
set chimerascr = ${PTSUITE}/chimera

if ($#argv < 1) then
   echo " "
   echo "script use:  RCSBvirusimages.csh  rcsbXXXXXX.cif rcsbXXXXXX.cif ..."
   echo "script generates both asymmetric unit and biological assembly images"
   echo "including structures split over multiple ids"
   echo " "
   exit 

else if (-f $1) then
    echo ''

else
    echo "Please check the file name. The name should be full rcsb id including the pdb extension"
    exit
endif

#get pdbid from from last field of HEADER and switch to lower case
# need to update this part -- get id from cif instead

set pdblist=''
set chimeraciflist=''
foreach i ($*)
  if (! -f $i) then
    echo $i 'file not found'
    exit
  endif
  set pdb = `awk '{ field = $NF }; NR==1 { print field }' $i | tr A-Z a-z`
  echo "pdb id read from file header (used for image file names): " $pdb
  set pdblist="$pdb $pdblist"
  set chimeraciflist="`pwd`/$i $chimeraciflist"
end

echo $chimeraciflist


#make biomt image
cp ${PTSUITE}/chimera/icos_cif.py .
echo " "
echo "making biological assembly image using BIOMT records"
echo " "
env CHIMERAPDB="$chimeraciflist"  chimera `pwd`/icos_cif.py  
echo " "

echo "making the biological assembly images for release"
foreach i ($*)
  set pdb = `awk '{ field = $NF }; NR==1 { print field }' $i | tr A-Z a-z`
  foreach size (500 250 80 65)
    convert -resize ${size}x${size} $i.jpg $pdb.pdb1-${size}.jpg 
  end
  display $pdb.pdb1-500.jpg &
end

rm -f icos_cif.py
rm -f icos_cif.pyc

#make crystal a.u. image
cp ${PTSUITE}/chimera/icos_cifncs.py .
echo " "
echo "making crystal asymmetric unit image using MTRIX records"
echo " "
env CHIMERAPDB="$chimeraciflist"  chimera ${PTSUITE}/chimera/icos_cifncs.py  
echo " "

echo "making the asymmetric unit images for release"
foreach i ($*)
  set pdb = `awk '{ field = $NF }; NR==1 { print field }' $i | tr A-Z a-z`
  foreach size (500 250 80 65)
  convert -resize ${size}x${size} $i.mtrix.jpg $pdb.pdb-${size}.jpg
  end
  display $pdb.pdb-500.jpg &
end


rm -f icos_cifncs.py
rm -f icos_cifncs.pyc

