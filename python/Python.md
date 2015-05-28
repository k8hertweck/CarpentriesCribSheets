##Software Carpentry Unix Instructor Crib Sheet
####Oklahoma State University, May 2015
####Kate L Hertweck, University of Texas at Tyler

####Before class:

* ipython notebook located in Dropbox folder, share link to .pynb with students so they can keep track
* enlarge text size
* check software installation: Python v 2.X (Anaconda)

####Checklist for class:
* class website: XXX
* etherpad: https://etherpad.mozilla.org/XXX
* use sticky notes

####SETUP

slides: http://swcarpentry.github.io/python-novice-inflammation/motivation.html 
why python?
We're teaching you how to program, and we have to use something
free, well documented, and everyone can run it
large userbase 
easy for novices to learn
data: CSV, comma separated values
goals: load data into memory, calculate average inflammation per day across all patients
lessons using Jupyter (ipython) notebooks
get data (if not already done)
http://swcarpentry.github.io/shell-novice/shell-novice-data.zip, move to Desktop, double click to unzip (if not already done), folder named “data”, can get there using: cd && cd Desktop/shell-novice/data
cd, git clone https://github.com/ouinformatics/osu-data 
orientation to notebook
terminal window, should start from directory with data and launch python, then ipython, then ipython notebook
create new notebook
go through orientation
saving and sharing with other people

####ANALYZING PATIENT DATA
Objectives: explain libraries, load library, read in data, assign values to variables, select data values, operations on arrays, simple graphs

libraries
import numpy
shift+enter to execute command
numpy.loadtxt(fname='inflammation-01.csv', delimiter=',')
first part is function call: function from particular library
dotted notation: numpy.loadtxt(...)
parameters: name of file, delimiter
assign variable: 
weight_kg = 55
print weight_kg
case sensitive
arithmetic: 
print 'weight in pounds:', 2.2 * weight_kg
separated with commas
change variable by assigning new
weight_kg =57.4
print 'weight in kilograms is now:', weight_kg
garbage 

####FILES AND DIRECTORIES
Objectives: paths, learn basic commands for working with files and directories, learn syntax of commands, tab-completion


####END CLASS
sticky note feedback
