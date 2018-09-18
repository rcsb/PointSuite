#!/bin/csh -f
# reset demo directories (remove output files)
if ($#argv < 1) then
   echo " "
   echo "script use:  cleandemo.csh DIRNAME or cleandemo.csh ALL"
   exit 
else if ($1 == "ALL"|| $1 == "all") then 
        set demo1 = "1CGM 1EI7 1IFD 1M4X 1RUG 2VF9 2W0C 2XD8 3N7X SPLIT IMPORT"
else if (-d $1) then
        set demo1 = $1
else
   echo "Please check directory name"
   exit
endif

foreach demodir ($demo1)
echo 'cleaning directory' ${demodir}
mkdir tmpdir
cp ${demodir}/${demodir}* tmpdir
rm -rf ${demodir}
mv tmpdir ${demodir}
end
