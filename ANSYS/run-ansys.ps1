# RUN ANSYS VIA POWERSHELL
# License (MIT):

# Copyright (c) 2016-2019 Christian Schubert

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#-------------------------------------------------------------------------------
##----------------Settings------------------------------------------------------

# Jobname
$jobname = 'apdl-example'

# Name of APDL script (must be located in working folder) to load
$skriptname = 'apdl_example.ans' 

# Working folder
$location = 'Z:\EXAMPLE\ANSYS'

# path for exchange (XC) folder
$xc_location = 'Z:\EXAMPLE\XC\'

# trials to run ansys (may be necessary if ansys behaves unstable...) 
$trys = 1

# Path to ansysxxx.exe
$ansysEXE = 'C:\Program Files\ANSYS Inc\v191\ansys\bin\winx64\ANSYS191'

# Ansys Settings

# Memory -m whole RAM ANSYS uses
# -db is the amount of RAM used for calculation if it is to small .page file is used which slows down calculations a lot ...
$memory = 20000
$dbmemory = 4000
$cpus = 4

# Parameter Liste
$paramList = '-j ' + $jobname + `
 ' -b ' + '-i ' + $skriptname +' -m ' + $memory + ' -db '+ $dbmemory `
+' -np ' + $cpus +' -noread ' + '-p aa_r ' + '-o ' + $jobname + '_OUT.TXT'


##----------------End Settings--------------------------------------------------


Set-Location -Path $location

echo 'Running Ansys Executable:' $ansysEXE 
echo 'With arguments: ' $paramList
echo ' '

for($i=1; $i -le $trys; $i++){
    echo 'starting Process...'
    Start-Process $ansysEXE -ArgumentList $paramList -WorkingDirectory $location -wait
    echo 'ended'
    echo ' '
}

echo 'Press [ENTER] to continue'
echo '...'
$userinput = [Console]::ReadLine()
