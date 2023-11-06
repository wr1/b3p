

REM This is a batch file to run the test case for the blade model 
REM http://www.calculixforwin.com/download/

set Y="blade_test.yml"
b3p --yml=%Y% clean build
b3p --yml=%Y% ccxprep --buckling=True
b3p --yml=%Y% ccxsolve forward -n 5  
b3p --yml=%Y% ccxpost
b3p --yml=%Y% ccxplot
