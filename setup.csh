#type "source setup.csh" to execute this script
set currentdir = `pwd` 
echo " "
echo "**** setting up pointsuite environment *******"
echo "the following lines are being executed:"
echo "(instead of running this script again"
echo " you can add them to your .(t)cshrc file):"
echo " "
echo "setenv PTSUITE" ${currentdir}
echo 'setenv PATH $PTSUITE/bin:$PATH'

setenv PTSUITE ${currentdir}
setenv PATH $PTSUITE/bin:$PATH
echo "*************** setup complete **************"

