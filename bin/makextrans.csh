#! /bin/csh -f
cat $1  > tmp
if ($2 == "ident") then
echo $3 ' "transform to crystal frame"  .' >> tmp
echo "1.0 0.0 0.0 0.0" >> tmp
echo "0.0 1.0 0.0 0.0" >> tmp
echo "0.0 0.0 1.0 0.0" >> tmp
else if (-e $2) then
  echo $3 ' "transform to crystal frame"  .' >> tmp
  set icount = 0
  foreach var (`cat $2`)
    if ($icount < 12) then
      echo  $var >> tmp
    endif
    @ icount = $icount + 1
  end
endif
mv tmp $1

