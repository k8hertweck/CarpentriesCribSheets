#!/usr/bin/env python3

#### Software Carpentry Programming with Python ####

#### Before class ####

# ipython notebook located in Dropbox folder, render in nbviewer, share link to latter so students can follow along (and download notebook)
# check software installation: Python v 3.X (Anaconda)

#### Objectives ####

# why python?
#   We're teaching you how to program, and we have to use something
#   free, well documented, and everyone can run it
#   large userbase
#   easy for novices to learn
#   super popular on campus!
# motivation: inflammation in patients who have been given new treatment for arthritis
#   load data into memory, calculate average inflammation per day across all patients, plot to share info with colleagues
#   data: CSV, comma separated values
#   rows contain information for a single patient (observations)
#   columns represent measurements on successive days

#### Setup ####

# many ways to interact with Python
#   python in terminal
#   ipython in terminal
#   save script in text editor
#   IDE like spyder
#   notebook: web application that combines code, graphs, and text
#   interactive mode in terminal, chevrons (>>>) is prompt, waiting for input
#   scripting mode: save commands in file (ends in .py), execute entire file at once
# about our tools
#   Anaconda: distribution (way of obtaining) Python;
#       includes extra packages like ipython, spyder
#   conda: package manager that comes with Anaconda, installs/updates packages
#   jupyter notebook: installed with Anaconda

# setting up jupyter project
#   launch Jupyter Notebook from Anaconda
#   terminal window must stay open, this is kernel (running python)
#   web browser is how you interact with notebook
#   create project directory (new folder), rename, then move into it
#   click "New" in upper right hand, then select "Python3"
#   creates notebook (*.ipynb, or ipython notebook file)
#   autosaves, or can save manually
#   click on title to rename
# executing code in a jupyter notebook:
#   enter code in cell and execute by pressing Shift + Return/enter
#   output is printed directly below cell, prefaced by Out[ ]:
#   add new cell with + button
#   can add Markdown cells with nicely formatted text
#   comments prefaced with # (not read/executed by python)
#   commands and output saved in notebook
#   talk about other menu options and buttons to remove/add/run cells
#   example notebook: https://github.com/rasilab/machkovech_2018/blob/master/scripts/NA43_competition.ipynb

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

# using variables
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

## Challenge: What values do the variables mass and age have after each statement in the following program? Test your answers by executing the commands.
mass = 47.5
age = 122
mass = mass * 2.0
age = age - 20
print(mass, age)

## Challenge: What does the following program print out?
first, second = 'Grace', 'Hopper'
third, fourth = second, first
print(third, fourth)

# libraries: collections of additional code that provide more functionality to perform specific tasks
# load library
import os
import urllib.request
import zipfile
import numpy

# download data
urllib.request.urlretrieve("http://swcarpentry.github.io/python-novice-inflammation/data/python-novice-inflammation-data.zip", "python-novice-inflammation-data.zip")
# unzip data
zipData = zipfile.ZipFile('python-novice-inflammation-data.zip')
zipData.extractall()

# load data into python (using library)
numpy.loadtxt(fname='data/inflammation-01.csv', delimiter=',')
# numpy.loadtxt(...) is a function call
#   run function loadtxt
#   belongs to numpy library
#   dotted notation in function call means whatever appears before dot contains the thing after the dot
# parameters: specific information that is sent to (passed) to function call
#   name of file
#   delimiter

# assign data to variable (so we can recall it later)
data = numpy.loadtxt(fname='data/inflammation-01.csv', delimiter=',')
# show the variable's value
print(data)
# what type of thing is data?
print(type(data))
# find type of data contained within array (data)
print(data.dtype)
# show shape of data
print(data.shape)
# output is rows, columns; rows are the individual patients, and the columns are their daily inflammation measurements
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

## Challenge: We can slice character strings as well! Given the following:
element = 'oxygen'
print('first three characters:', element[0:3])
print('last three characters:', element[3:6])
# What is the value of element[:4]? What about element[4:]? Or element[:]?
element[4:]
element[:]

## Challenge: What is element[-1]? What is element[-2]?
element[-1]
element[-2]

# Given those answers, explain what element[1:-1] does.
# Creates a substring from index 1 up to (not including) the final index, effectively removing the first and last letters from ‘oxygen’

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

# view max inflammation per patient or per day
# create temporary array for data desired
patient_0 = data[0, :] # 0 on the first axis (rows), everything on the second (columns)
print('maximum inflammation for patient 0:', numpy.max(patient_0))

# combine selection and function call (skip temp variable)
print('maximum inflammation for patient 2:', numpy.max(data[2, :]))

# average across all rows (axis 0)
print(numpy.mean(data, axis=0))
# confirm shape of array
print(numpy.mean(data, axis=0).shape)
# average across all columns (axis 1): avg inflammation per day for all patients
print(numpy.mean(data, axis=1))

## Challenge: 

#### Visualizing data ####

# make pylot available from matplotlib (de facto plotting library)
import matplotlib.pyplot
# allow plots to appear when using show()
%matplotlib inline

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

## Challenge: Create a plot showing the standard deviation (numpy.std) of the inflammation data for each day across all patients.
std_plot = matplotlib.pyplot.plot(numpy.std(data, axis=0))
matplotlib.pyplot.show()

# grouping plots: complete set of code

# load libraries
import numpy
import matplotlib.pyplot

# load data from file
data = numpy.loadtxt(fname='data/inflammation-01.csv', delimiter=',')

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

## Challenge: Modify the program to display the three plots on top of one another instead of side by side.
import numpy
import matplotlib.pyplot

data = numpy.loadtxt(fname='inflammation-01.csv', delimiter=',')

# change figsize (swap width and height)
fig = matplotlib.pyplot.figure(figsize=(3.0, 10.0))

# change add_subplot (swap first two parameters)
axes1 = fig.add_subplot(3, 1, 1)
axes2 = fig.add_subplot(3, 1, 2)
axes3 = fig.add_subplot(3, 1, 3)

axes1.set_ylabel('average')
axes1.plot(numpy.mean(data, axis=0))

axes2.set_ylabel('max')
axes2.plot(numpy.max(data, axis=0))

axes3.set_ylabel('min')
axes3.plot(numpy.min(data, axis=0))

fig.tight_layout()

matplotlib.pyplot.show()

## Challenge: How would you alter the limits on the x and y axes?
axes3.set_ylim(0,6)
# A more automated approach
min_data = numpy.min(data, axis=0)
axes3.set_ylabel('min')
axes3.plot(min_data)
axes3.set_ylim(numpy.min(min_data), numpy.max(min_data) * 1.1)

#### Repeating actions with loops ####

# Objectives: write for loop to repeat simple actions, trace changes to variables

# what if we wanted to print each character in a word on a line of its own?
word = 'lead'
print(word[0])
print(word[1])
print(word[2])
print(word[3])
# try a different word
# this doesn't scale well, and is fragile (creates error if word is shorter, doesn't print all for longer word)

# print with for loop
for char in word:
    print(char)
# syntax:
#   for variable in collection:
#       do things using variable
# note on choosing meaningful variable names

# for loop that repeatedly updates variable
length = 0 # define external variable
for vowel in 'aeiou': # initialize for loop for string (not variable)
    length = length + 1 # define continuously updated variable internal to loop
print('There are', length, 'vowels') # report output at end of loop

# loop variables still exist after loop ends!
letter = 'z'
for letter in 'abc':
    print(letter)
print('after the loop, letter is', letter)

# finding length of a string is a built-in function!
print(len('aeiou'))

## Challenge:

#### Storing multiple values in lists ####

# Objectives: create and index lists of simple values, change values of elements, append to list, reorder and slice lists, create and manipulate nested lists

# create a list
odds = [1, 3, 5, 7]
# recall list
print('odds are:', odds)
# select individual elements via indexing
print('first and last:', odds[0], odds[-1])
# loop over list
for number in odds:
    print(number)
# you can change values in a list (but not individual characters in a string)
# list example
names = ['Curie', 'Darwing', 'Turing']  # typo in Darwin's name
print('names is originally:', names)
names[1] = 'Darwin'  # correct the name
print('final value of names:', names)
# string example
name = 'Darwin'
name[0] = 'd'

# two variables can refer to the same list; modifying one modifies both!
salsa = ['peppers', 'onions', 'cilantro', 'tomatoes']
my_salsa = salsa        # <-- my_salsa and salsa point to the *same* list data in memory
salsa[0] = 'hot peppers'
print('Ingredients in my salsa:', my_salsa)

# a better way is to make a copy of the original list
salsa = ['peppers', 'onions', 'cilantro', 'tomatoes']
my_salsa = list(salsa)        # <-- makes a *copy* of the list
salsa[0] = 'hot peppers'
print('Ingredients in my salsa:', my_salsa)

# lists can contain other lists
x = [['pepper', 'zucchini', 'onion'],
     ['cabbage', 'lettuce', 'garlic'],
     ['apple', 'pear', 'banana']]
# print first row
print([x[0]])
print(x[0])
# print first item in first row
print(x[0][0])

# lists can contain elements of different types
sample_ages = [10, 12.5, 'Unknown']

# append to the list
odds.append(11)
print('odds after adding a value:', odds)

# remove first element
del odds[0]
print('odds after removing the first element:', odds)

# reverse the list
odds.reverse()
print('odds after reversing:', odds)

# slicing lists
binomial_name = "Drosophila melanogaster"
group = binomial_name[0:10]
print("group:", group)

species = binomial_name[11:24]
print("species:", species)

chromosomes = ["X", "Y", "2", "3", "4"]
autosomes = chromosomes[2:5]
print("autosomes:", autosomes)

last = chromosomes[-1]
print("last:", last)

# omit the first or last range to indicate the start or end
date = "Monday 4 January 2016"
day = date[0:6]
print("Using 0 to begin range:", day)
day = date[:6]
print("Omitting beginning index:", day)

months = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
sond = months[8:12]
print("With known last position:", sond)
sond = months[8:len(months)]
print("Using len() to get last entry:", sond)
sond = months[8:]
print("Omitting ending index:", sond)

#### Analyzing data from multiple files ####

# import library to find files/directories
import glob

# get names of all csv files in current directory
print(glob.glob('inflammation*.csv'))

# combine previous content together and analyze all files

import numpy
import matplotlib.pyplot

filenames = sorted(glob.glob('inflammation*.csv'))
filenames = filenames[0:3]
for f in filenames:
    print(f)

    data = numpy.loadtxt(fname=f, delimiter=',')

    fig = matplotlib.pyplot.figure(figsize=(10.0, 3.0))

    axes1 = fig.add_subplot(1, 3, 1)
    axes2 = fig.add_subplot(1, 3, 2)
    axes3 = fig.add_subplot(1, 3, 3)

    axes1.set_ylabel('average')
    axes1.plot(numpy.mean(data, axis=0))

    axes2.set_ylabel('max')
    axes2.plot(numpy.max(data, axis=0))

    axes3.set_ylabel('min')
    axes3.plot(numpy.min(data, axis=0))

    fig.tight_layout()
    matplotlib.pyplot.show()

## Challenge:

#### Making choices ####

# Objectives: write conditional statements including if, elif, else; evaluate expressions containing and and or

# tell python to take different actions with if statements
num = 37
if num > 100:
    print('greater')
else:
    print('not greater')
print('done')

# don't need else; can also do nothing
num = 53
print('before conditional...')
if num > 100:
    print(num,' is greater than 100')
print('...after conditional')

# can have multiple alternatives using elif
num = -3

if num > 0:
    print(num, 'is positive')
elif num == 0: # double equal sign is necessary; single used to assign values
    print(num, 'is zero')
else:
    print(num, 'is negative')

# combine tests using and, when both parts must be true
if (1 > 0) and (-1 > 0):
    print('both parts are true')
else:
    print('at least one part is false')

# combine tests using or if at least one part must be true
if (1 < 0) or (-1 < 0):
    print('at least one test is true')

# checking for problems in inflammation data
import numpy # if not already done

# check if max inflammation equals day number (error in data entry)
max_inflammation_0 = numpy.max(data, axis=0)[0]
max_inflammation_20 = numpy.max(data, axis=0)[20]

if max_inflammation_0 == 0 and max_inflammation_20 == 20:
    print('Suspicious looking maxima!')

# check if mins are all zero (healthy patient)
if numpy.sum(numpy.min(data, axis=0)) == 0:
    print('Minima add up to zero!')

# combine together with data
data = numpy.loadtxt(fname='inflammation-01.csv', delimiter=',')

max_inflammation_0 = numpy.max(data, axis=0)[0]
max_inflammation_20 = numpy.max(data, axis=0)[20]

if max_inflammation_0 == 0 and max_inflammation_20 == 20:
    print('Suspicious looking maxima!')
elif numpy.sum(numpy.min(data, axis=0)) == 0:
    print('Minima add up to zero!')
else:
    print('Seems OK!')

# test on another dataset

## Challenge:
