##Software Carpentry Unix Instructor Crib Sheet
####Kate L Hertweck, University of Texas at Tyler

####Before class:
* set up RStudio
* set up Socrative questions 
	* teacher room: PAW5AYWM
  	* start quiz, teacher paced, disable student names
* put history in accessible place for students

**Intro to Software Carpentry:** http://swcarpentry.github.io/slideshows/introducing-software-carpentry/index.html#slide-0 

####Checklist for class:
* class website: http://XXX
* etherpad: https://etherpad.mozilla.org/XXX
* Kate's history: https://www.dropbox.com/XXX
* instructor's guide for R: http://swcarpentry.github.io/r-novice-inflammation/instructors.html
* have zip of data linked somewhere for students to download

####SETUP (5 min)

* we're going to learn R for data analysis!
* real goal not to teach R, but teach basic concepts that you can use in any language
* why R?
	* we have to use something
	* it's free, well-documented, and runs on most machines
	* active community, especially among scientists
	* lots of code already written (and available for reuse!) for many different types of analyses
* download data and put in SWC workshop folder: http://swcarpentry.github.io/python-novice-inflammation/python-novice-inflammation-data.zip
	* don't worry if you see "python" in name, unzip and it should be called "data"
	* can embed link in workshop webpage for students to download

####INTRO TO RSTUDIO (10 min)

* RStudio is just a handy interface to work with R (some folks may have used the regular R console already)
	* R works the same, regardless of whether you're using the R console, R GUI, or RStudio
* running commands in console
	* when R is ready to accept a command, you'll see `>` as the command prompt in the console
	* you can perform math: `5 + 5`
	* hit `enter` to execute commands
	* answer appears in console
	* if you see a `+` instead of the command prompt, R is expecting you to finish the previous command. This is probably an error caused by missing parentheses, quotation marks, etc. You can try entering one of those symbols until an error message appears, or hit `esc` (although sometimes that makes RStudio crash)
* get project set up
	* File -> New Project -> Existing Directory (where data are stored), creates a few files to that will make managing your work easier, some windows may change
	* File -> New File -> R script (save as something like `workshop.R`, `.R` lets you know it's R code)
* parts of RStudio
	* top left: R script, where you can type and save commands for later reference (or create executable scripts)
	* bottom left: console, where commands are executed and you can see resulting output
	* top right: environment, which lets you know what pieces of information R "remembers"
	* bottom right: lots of things! file hierarchy, plots, packages, help
	* everything can be resized or minimized

####PROGRAMMING WITH R
**Objectives:** read in tabular data, assign variables, select individual values/subsections, perform operations on data frames, display simple graphs

* csv: rows and columns, separated by commas
* exists as an independent file. we want to:
	* load in data
	* calculate average inflammation per day across all patients
	* plot results
* loading data
	* `setwd("~/swc)` tell computer where files are
	* enter this in the R script file if you want to be able to save it for later
	* to execute from R script, hold `control` and hit `enter`
	* can also use menu option, Session -> Set Working Directory -> Choose Directory (but this doesn't get saved in script!)
	* once your working directory matches where the file is located, you can enter data
	* `read.csv(file="data/inflammation-01.csv, header = FALSE)`
	* `read.csv(...)` is a function call that asks R to run function `read.csv` with certain parameters
	* arguments are the name of the file we want to read and whether the first line of the file contains a header row
	* there are lots of arguments for this function (can find using `?read.csv` or using the help tab in the bottom right window
	* the result of this command was to print the file to the console. This doesn't mean R "remembers" it!
	* assign data to variable to retain it in R for future use
	* `weight_kg <- 55` then `control + enter`
	* spaces are not strictly necessary, but good programming practice for readability
	* now when you use `weight_kg`, R will recall the value you assigned
	* do arithmetic: `2.2 * weight_kg`
	* add comment above this in R script saying what it does: `# weight in pounds`
	* save as variable: `weight_lb <- 2.2 * weight_kg`, then `weight_lb`
	* change object's value: `weight_kg <- 57.5`
	* what do you expect if you type `weight_kg`?
	* if you reassign the value, R will "forget" the old value
	* now we're ready to save data to file! 
	* `dat <- read.csv(file="data/inflammation-01.csv, header = FALSE)`
	* no output, but you should see something appear in environment in upper right window
	* try `head(dat)`: prints first few linse of data
	* **Socrative challenge question:** What is the output from the following code:
	```
	mass <- 47.5
	age <- 122
	mass <- mass * 2.0
	age <- age - 20
	```
* manipulating data
	* what type of thing is `dat`? `class(dat)`
	* data frame: like a spreadsheet from MS Excel, common for experimental data with individual observations in rows and variables in columms
	* what are the dimensions (shape)? `dim(dat)`
	* always references rows, columns
	* extract particular values:
	```
	# first value in dat
	dat[1, 1]
	# middle value in dat
	dat[30, 20]
	```
	* select sections of data using ranges of rows and columns (slice)
	```
	# select section of data
	dat[1:4, 5:10]
	```
	* select non-contiguous values using combine command:
	```
	# select non-contiguous values
	dat[c(3, 8, 37, 56), c(10, 14, 29)]
	```
	* leave out rows or columns
	``` 
	# select all columns but certain rows
	dat[5, ]
	# select all rows but certain columns
	dat[, 16]
	```
	* perform manipulations on data
	```
	# first row, all columns
	patient_1 <- dat[1, ]
	# max inflammation for patient 1
	max(patient_1)
	```
	* combine data selection and function calls:
	``` 
	# max inflammation for patient 2
	max(dat[2, ])
	```
	* other funtions:
	```
	# minimum inflammation on day 7
	min(dat[, 7])
	# mean inflammation on day 7
	mean(dat[, 7])
	# median inflammation on day 7
	median(dat[, 7])
	# standard deviation of inflammation on day 7
	sd(dat[, 7])
	```
	* what about finding maximum inflammation for all patients? or original task: average for each day? use `apply`, which allows repeat of function on rows (1) or columns (2), called the margin
	```
	# repeat function across rows
	avg_patient_inflammation <- apply(dat, 1, mean)
	# repeat function across columns
	avg_day_inflammation <- apply(dat, 2, mean)
	```
	* can also arrow up in console to recall previous commands, but remember that it won't be saved in your R script
	* some of these common operations also have their own functions (i.e., `rowMeans` and `colMeans` (remember that capitalization matters!)
	* **Socrative challenge question** (subsetting data)
* plotting data
	* basic intro to R plotting features
	* `plot(avg_day_inflammation)` to give R a vector of data with numbers corresponding to average inflammation per day
	* y-axis is average inflammation level, x-axis is order (index) of values in vector
	* this pattern is a bit suspicious! we would expect a sharper rise and slower fall
	* try other statistics
	```
	# maximum inflammation per day
	max_day_inflammation <- apply(dat, 2, max)
	plot(max_day_inflammation)
	# minimum inflammation per day
	min_day_inflammation <- apply(dat, 2, min)
	plot(min_day_inflammation)
	```
	* this lets us diagnose a problem with the data
	* **Socrative challenge question**: create plot showing standard deviation for each day across all patients

####END CLASS
* remove R history from dropbox (remind students to save first if they want it)
* give students links to complete R lessons (and intermediate lessons)
* sticky note feedback
