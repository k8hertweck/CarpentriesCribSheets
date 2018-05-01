## Data Carpentry R Ecology Instructor Crib Sheet

### Before class (instructor):
* set up RStudio
* set up Socrative questions 
	* teacher room: PAW5AYWM
  	* start quiz, teacher paced, disable student names
* student resources:
	* Instructor history in accessible place for students (e.g., Dropbox); share link
	* class website
	* etherpad
	* check installed software: R and RStudio

## BEFORE WE START
* **Objectives**: Learn why R, intro to RStudio, file hierarchy, R syntax, get help

### Basics of R

* R is a useful tool for both statistics and data science
* why R?
	* we have to use something
	* it's free, well-documented, and runs on most machines
	* active community, especially among scientists
	* lots of code already written (and available for reuse!) for many different types of analyses

### Intro to R and RStudio

* RStudio is just a handy interface to work with R (some folks may have used the regular R console already)
	* R works the same, regardless of whether you're using the R console, R GUI, or RStudio
* running commands in console
	* when R is ready to accept a command, you'll see `>` as the command prompt in the console
	* you can perform math: `5 + 5`
	* hit `enter` to execute commands
	* answer appears in console
	* if you see a `+` instead of the command prompt, R is expecting you to finish the previous command. This is probably an error caused by missing parentheses, quotation marks, etc. You can try entering one of those symbols until an error message appears, or hit `esc` (although sometimes that makes RStudio crash)
	
### Setting up a project

* get project set up
	* File -> New Project -> New directory -> Empty project, name `data-carpentry`, creates a few files to that will make managing your work easier, some windows may change
	* File -> New File -> R script (save as something like `data-carpentry-script.R`, `.R` lets you know it's R code)
* parts of RStudio
	* top left: R script, where you can type and save commands for later reference (or create executable scripts)
	* bottom left: console, where commands are executed and you can see resulting output
	* top right: environment, which lets you know what pieces of information R "remembers"
	* bottom right: lots of things! file hierarchy, plots, packages, help
	* everything can be resized or minimized

### Organizing working directory
* original data in `data/`
* intermediate datasets in `data_output/`
* figures in `figure_output/`

### Getting help
* `??kruskal` to search throughout help for all functions
* `sessionInfo()` to find information about how R is set up for you.
* online resources

## INTRO TO R
* **Objectives**: expand knowledge of R syntax, objects/assignments, vector/data types

### Creating objects in R
* start modeling saving commands in R script and executing with keyboard shortcut into console
* `#` is a comment in R; anything to the right on line not executed by console (just printed out)
* link to other keyboard shortcuts 
* type math in console: `3 + 5`
* this prints the result to screen, not "remembered" by R
* `weight_kg <- 55`
* `<-` is assignment operator
* this variable (object) shows up in the environment window
* object names:
	* should be descriptive, but not too long
	* can't start with number
	* should be different than other pre-defined commands (mean, data, etc)
	* avoid dots
* R doesn't print anything to console when creating variable
* can recall using `weight_kg` or `(weight_kg <- 55)`
* now object is in memory, we can perform other operations on it
* `2.2 * weight_kg`
* assign a new value to variable: `weight_kg <- 57.5`
* `2.2 * weight_kg` gives updated value
* `weight_lb <- 2.2 * weight_kg` create new object 
* `weight_kg <- 100` send new value to object 
* what is `weight_lb` now?
* **Challenge question (Socrative)** What are the values after each statement in the following?
```
mass <- 47.5            # mass?
age  <- 122             # age?
mass <- mass * 2.0      # mass?
age  <- age - 20        # age?
mass_index <- mass/age  # mass_index?
```
* Answer: mass 95, age 102

### Functions and arguments
* functions are little programs in R that let you perform helpful tasks
* `round(3.14159)`
* learn what other options are available for this function: `args(round)`
* fund more information: `?round`
* apply this argument for this function: `round(3.14159, digits=2)`
* names and order can be important: if not named, need to be in right order, if named, can put in any order: `round(digits=2, x=3.14159)`


### Vectors and data types
* vector: basic data structure in R
* groups of values (numbers or characters)
* assign in list: `weights <- c(50, 60, 65)`
* recall with `weights`
* can also be characters: `animals <- c(“mouse”, “rat”, “dog”, "cat")`
* inspecting contents of vector:
	* how many elements? `length(weights)` and `length(animals)`
	* type of element? `class(weights)` and `class(animals)`
	* overview? `str(weights)` and `str(animals)`
* add elements using `c()` function (like concatenate)
	* `weights <- c(weights, 90) # add at end`
	* `weights <- c(30, weights) # add at beginning`
	* this may be useful for building a dataset, or autoupdating as we program
* data types
	* numeric
	* character
	* logical (true/false)
	* integer
	* complex

### Subsetting vectors 
* extract individual values from within vector: `animals[2]`
* extract multiple values: `weights[c(1,3)]`
* other types of data structures:
	* list
	* matrix
	* data.frames
	* factors
* **Challenge questions** 
	* We’ve seen that atomic vectors can be of type character,
  numeric, integer, and logical. But what happens if we try to mix these types in
  a single vector?
		* _Answer_: R implicitly converts them to all be the same type
	* **Question**: Why do you think it happens?
		* _Answer_: Vectors can be of only one data type. R tries to convert (=coerce)
  the content of this vector to find a "common denominator".
	* **Question**: Can you draw a diagram that represents the hierarchy of the data
  types?
		* _Answer_: `logical --> numeric --> character <-- logical`

### Conditional subsetting
* `weight_g > 50    # will return logicals with TRUE for the indices that meet the condition`
* `weight_g[weight_g > 50] # select only values above 50`
* multiple tests: `weight_g[weight_g < 30 | weight_g > 50] # OR, can also use & for AND`

### Missing data
* missing data are represented as NA in R
* use `na.rm=TRUE` to ignore missing data (otherwise might get NA as answer)
```
heights <- c(2, 4, 4, NA, 6)
mean(heights)
max(heights)
mean(heights, na.rm = TRUE)
max(heights, na.rm = TRUE)
```
* also `is.na()`, `na.omit()`, and `complete.cases()`
