#!/usr/bin/env python3

#### Software Carpentry Programming with Python ####

#### Before class ####

# ipython notebook located in Dropbox folder, share link to .pynb with students so they can keep track
# enlarge text size
# check software installation: Python v 3.X (Anaconda)

#### Objectives ####

# slides: http://swcarpentry.github.io/python-novice-inflammation/motivation.html
# why python?
#   We're teaching you how to program, and we have to use something
#   free, well documented, and everyone can run it
#   large userbase
#   easy for novices to learn
# data: CSV, comma separated values
# motivation: load data into memory, calculate average inflammation per day across all patients

#### Setup ####

# lessons using Jupyter (ipython) notebooks
# get data (if not already done):
#   http://swcarpentry.github.io/shell-novice/shell-novice-data.zip, move to Desktop, double click to unzip (if not already done), folder named “data”, can get there using: cd && cd Desktop/shell-novice/data
#   cd, git clone https://github.com/ouinformatics/osu-data
# orientation to notebook
#   terminal window, should start from directory with data and launch python, then ipython, then ipython notebook
# create new notebook
# go through orientation
# saving and sharing with other people

#### Analyzing patient data ####

# Objectives: intro to libraries, read in data, assign values to variables, select data values, operations on arrays, simple graphs

# include human-readable (but not python interpreted) comments following hash signs

# use python as calculator
3 + 5
# shift+enter to execute command

# assign value to variable
weight_kg = 60
# variable names:
#   can include letters, digits, and underscores
#   cannot start with a digit
#   are case sensitive

# data types:
#   integers
#   floating point numbers (decimals)
#   strings (characters)
# weight_kg is integer
# to create as floating point
weight_kg = 60.0
# to create string
weight_kg_text = 'weight in kilograms:'

# display value of variable
print(weight_kg)
# print is a function
# can display multiple items at once
print(weight_kg_text, weight_kg)

# perform arithmetic inside print function
print('weight in pounds:', 2.2 * weight_kg)
# note: this doesn't change the value of weight_kg!
print(weight_kg)

# assign new value to weight_kg
weight_kg = 65.0
print('weight in kilograms is now:', weight_kg)
# variable names as sticky notes (analogy)

## Challenge:

# libraries: describe what they are
# load library
import numpy

# load data into python (using library)
numpy.loadtxt(fname='inflammation-01.csv', delimiter=',')
# numpy.loadtxt(...) is a function call
#   run function loadtxt
#   belongs to numpy library
#   dotted notation in function call means whatever appears before dot contains the thing after the dot
# parameters: specific information that is sent to (passed) to function call
#   name of file
#   delimiter

# assign data to variable (so we can recall it later)
data = numpy.loadtxt(fname='inflammation-01.csv', delimiter=',')
# show the variable's value
print(data)
# what type of thing is data?
print(type(data))
# find type of data contained within array (data)
print(data.dtype)
# show shape of data
print(data.shape)
# output is rows, columns
# arrays have members, or attributes, which use the dot nomenclature because they have the same part-and-whole relationship

# access a specific value
print('first value in data:', data[0, 0])
# python begins indexing (counting) at 0
print('middle value in data:', data[30, 20])

# select sections of data (slicing)
print(data[0:4, 0:10])
# end bound is NOT inclusive (up to but not including)
# can start at indeces besides 0
print(data[5:10, 0:10])
# use empty bound to include the end of axis
small = data[:3, 36:] # assign to value
print('small is:')
print(small)

# perform math on array
doubledata = data * 2.0
# view output
print('original:')
print(data[:3, 36:])
print('doubledata:')
print(doubledata[:3, 36:])

# perform operation involving two arrays
tripledata = doubledata + data
print('tripledata:')
print(tripledata[:3, 36:])

## Challenge:

# perform calculation across entire array
print(numpy.mean(data)) # find mean

# use multiple assignment to obtain descriptive values of data
maxval, minval, stdval = numpy.max(data), numpy.min(data), numpy.std(data)

print('maximum inflammation:', maxval)
print('minimum inflammation:', minval)
print('standard deviation:', stdval)

# to find available functions and information:
#   type name of something, followed by dot, then hit tab (ipython and notebooks)
#   select a function or attribute and add question mark to find help documentation
#   help(thing.attribute) is same as above

# create temporary array for data desired
patient_0 = data[0, :] # 0 on the first axis (rows), everything on the second (columns)
print('maximum inflammation for patient 0:', numpy.max(patient_0))

# combine selection and function call (skip temp variable)
print('maximum inflammation for patient 2:', numpy.max(data[2, :]))

# average across all rows (axis 0)
print(numpy.mean(data, axis=0))
# confirm shape of array
print(numpy.mean(data, axis=0).shape)
# average across all columns (axis 1)
print(numpy.mean(data, axis=1))

## Challenge:

#### Visualizing data ####

# make pylot available from matplotlib (de facto plotting library)
import matplotlib.pyplot
# allow plots to appear when using show()
%matplotlib inline # % only applies to functions valid in notebook environment

# make heatmap from data
image = matplotlib.pyplot.imshow(data)
matplotlib.pyplot.show()

# plot average inflammation over time
ave_inflammation = numpy.mean(data, axis=0)
ave_plot = matplotlib.pyplot.plot(ave_inflammation)
matplotlib.pyplot.show()

# plot max over time
max_plot = matplotlib.pyplot.plot(numpy.max(data, axis=0))
matplotlib.pyplot.show()

# plot min over time
min_plot = matplotlib.pyplot.plot(numpy.min(data, axis=0))
matplotlib.pyplot.show()

# grouping plots: complete set of code

# load libraries
import numpy
import matplotlib.pyplot

# load data from file
data = numpy.loadtxt(fname='inflammation-01.csv', delimiter=',')

# create space to place plot
fig = matplotlib.pyplot.figure(figsize=(10.0, 3.0)) # state dimensions of figure

# add subplots, parameters are number of subplots, number of columns, which subplot (left-to-right, top-to-bottom)
axes1 = fig.add_subplot(1, 3, 1)
axes2 = fig.add_subplot(1, 3, 2)
axes3 = fig.add_subplot(1, 3, 3)

# title axes
axes1.set_ylabel('average')
axes1.plot(numpy.mean(data, axis=0))

axes2.set_ylabel('max')
axes2.plot(numpy.max(data, axis=0))

axes3.set_ylabel('min')
axes3.plot(numpy.min(data, axis=0))

# spread out graphs
fig.tight_layout()

# show plot
matplotlib.pyplot.show()

## Challenge:
