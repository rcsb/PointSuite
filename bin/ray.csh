#/bin/csh -f
# get rid of specific cif categories in a file using command-line tools (grep, sed, awk)
# for loops, file is processed in two parts and then put back together

set filein = ${1}

# need to use quotations and append backslash for inclusion of period in grep
foreach ciftag ("_pdbx_point_symmetry\." "_pdbx_helical_symmetry\." "_pdbx_struct_assembly\." "_pdbx_struct_assembly_gen\." "_pdbx_struct_oper_list\.")
#foreach ciftag ( "_pdbx_struct_assembly\." )

#find out if category is present
if (`grep -c $ciftag $filein` == 0) then
   echo 'no' $ciftag 'category present'
else
#find out if it is a loop
grep -B 1 -m 1 $ciftag $filein >tmp 

if (`grep -c "loop_" tmp` == 0) then
   #not a loop
   echo 'no loop detected for' $ciftag
   grep -v $ciftag $filein > tmp
   mv tmp edited_file.cif
   set filein = 'edited_file.cif'

else
   #loop
   echo 'loop detected for' $ciftag
   #top file part
   sed -n '1,/'$ciftag'/ p' $filein | sed '$ d' | sed '$ d'  > start
   # bottom part starts at next loop or cif tag
   awk '/'$ciftag'/,0' $filein | grep -v $ciftag |  awk '/_[a-zA-Z]\+\.[a-zA-Z]|loop_/,0' > end
   cat start end > tmp
   mv tmp edited_file.cif
   set filein = 'edited_file.cif'
endif
endif
end


