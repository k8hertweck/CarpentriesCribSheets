## Data Carpentry R Instructor Crib Sheet
### Kate L Hertweck, University of Texas at Tyler

### Before class (instructor):
* set up RStudio
* set up Socrative questions 
	* teacher room: PAW5AYWM
  	* start quiz, teacher paced, disable student names
* put history in accessible place for students (e.g., Dropbox)

### Checklist for class (student resources):
* class website: XXX
* etherpad: XXX
* Kate's history: XXX
* installed software: R and RStudio

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

### Functions and arguments
* start modeling saving commands in R script and executing with keyboard shortcut into console
* `#` is a comment in R; anything to the right on line not executed by console (just printed out)
* link to other keyboard shortcuts 
* functions are little programs in R that let you perform helpful tasks
* `round(3.14159)`
* learn what other options are available for this function: `args(round)`
* fund more information: `?round`
* apply this argument for this function: `round(3.14159, digits=2)`
* names and order can be important: if not named, need to be in right order, if named, can put in any order: `round(digits=2, x=3.14159)`

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

### The R syntax
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

## STARTING WITH DATA

### Survey data
* .csv file containing species and weight
	* each row has info for one collection
	* each column contains info collected for each specimen
* download data and save in data directory: `download.file("https://ndownloader.figshare.com/files/2292169", "data/portal_data_joined.csv")`
* load data: `surveys <- read.csv('data/portal_data_joined.csv')`
	* no output, but data are saved in environment
* recall with `surveys`
* limited output: `head(surveys)`
* data.frame: table where columns are vectors of same length
* `str(surveys)`
* **Challenge question** Based on the output of `str(surveys)`:
	* What is the class of the object `surveys`?
	* How many rows and how many columns are in this object?
	* How many species have been recorded during these surveys?

### Factors
* factors represent categorical data, can be ordered or unordered
* contain pre-defined set of values known as levels (assigned alphabetically)
* `sex <- factor(c("male", "female", "female", "male"))`
* `levels(sex)` and `nlevels(sex)`
* may need to specify the levels
* `food <- factor(c("low", "high", "medium", "high", "low", "medium", "high"))`
* `levels(food)`
* `food <- factor(food, levels=c("low", "medium", "high"))`
* `levels(food)`
* `min(food) ## doesn't work`
* `food <- factor(food, levels=c("low", "medium", "high"), ordered=TRUE)` 
* `levels(food)`
* `min(food) ## works!`

### Converting factors
* converting factors to characters only requires `as.character()`
* converting factors to numeric is tricker: `f <- factor(c(1, 5, 10, 2))`
* `as.numeric(f) ## wrong! and there is no warning...`
* `as.numeric(as.character(f)) ## works...`
* `as.numeric(levels(f))[f] ## The recommended way.`
	* obtain all the factor levels using `levels(f)`
	* convert these levels to numeric values using `as.numeric(levels(f))`
	* access these numeric values using the underlying integers using `[f]`

## THE DATA.FRAME CLASS
* **Objectives**: understand data.frames, sequences, access elements of data.frame

### What are data.frames?
* default data structure for tabular data, used for statistics and plotting
* import using `read.csv()` or `read.table()`; many default options
* can force words to be characters (instead of factors) by setting `stringsAsFactors=FALSE`
* import to know how R is interpreting your data!
* **Challenge question** There are mistakes in this data.frame. How to fix?
```
author_book <- data.frame(author_first=c("Charles", "Ernst", "Theodosius"),
                          author_last=c(Darwin, Mayr, Dobzhansky),
                          year=c(1942, 1970))
```
* Answer: author_last in quotes, add NA for first cell in year

### Inspecting data frames
* size:
	* `dim()` number of rows, number of columns
	* `nrow()` rows
	* `ncol()` columns
* content
	* `head()`
	* `tail()`
* names
	* `names()` column names
	* `rownames()` row names
* summary
	* `str()` structure, class, length, content (by column)
	* `summary()` summary stats for each column

### Indexing and sequences
* you can extract values from a vector using brackets
* recall animals `animals <- c("mouse", "rat", "dog", "cat")`
* `animals[2]`
* `animals[c(3,2)]`
* `animals[2:4] #colon selects all items between`
* `more_animals <- animals[c(1:3,2:4)]
* R starts counting at one (this is not true of all languages!), called indexing
* select based on more complex patterns: 
	* `seq(1, 10, by=2) #counts by two within range, starting with first value`
	* `seq(5, 10, length.out=3) #splits interval into three equal values`
	* `seq(50, by=5, length.out=10) #counts by 5 starting at 50, returns 10 values`
	* `seq(1, 8, by=3) # sequence stops to stay below upper limit`
* for data frames, you may need to reference both rows and columns
	* R always assumes first value is row, second is column
	* `surveys[1] # first column`
	* `surveys[1, 1] # first row, first column`
	* `surveys[1, 6] # first row, sixth column`
	* `surveys[1:3, 7] # first three elements`
	* `surveys[3, ] # the 3rd element for all columns`
	* `surveys[, 8] # the entire 8th column`
	* `head_surveys <- surveys[1:6, ] # surveys[1:6, ] is equivalent to head(surveys)`
* can also call columns by name
	* `surveys[, "species_id"]`
	* `surveys[["species_id"]]`
	* `surveys$species_id`
	* for our purposes, these are equivalent
	* the last also includes partial matching, so `surveys$d` includes `day` column
	* RStudio has autocomplete
* **Challenge question**

### Conditional subsetting
* logical vectors: assigning value to pieces of data
```
animals <- c("mouse", "rat", "dog", "cat")
animals[c(TRUE, FALSE, TRUE, TRUE)]
```
* test the logical assignments:
```
animals != "rat"
animals[animals != "rat"]
animals[animals == "cat"]
```
* combine multiple tests using `&` (both conditions are true, AND) or `|`
(at least one of the conditions if true, OR):
```
animals[animals == "cat" & animals == "rat"] # returns nothing
animals[animals == "cat" | animals == "rat"] # returns both rat and cat
```
* function `%in%` allows you to test if a value if found in a vector
```
animals %in% c("rat", "cat", "dog", "duck")
animals[animals %in% c("rat", "cat", "dog", "duck")]
```
* test whether the elements of your vector are less than or greater than a given value:
```
dates <- c(1960, 1963, 1974, 2015, 2016)
dates >= 1974
dates[dates >= 1974]
dates[dates > 1970 & dates <= 2015]
dates[dates < 1975 | dates > 2016]
```
* **Challenge** Why does `"four" > "five"` returns `TRUE`?
	* Answer: this is based on alphabet

## AGGREGATING AND ANALYZING DATA WITH DPLYR
* **Objectives**: describe packages, selecting columns and filtering rows, pipes, mutate, split/apply/combine

### Data manipulation with dplyr
* packages: collections of additional functions to help you perform more operations
* while some commands are included with a general R installation ("base R"), individual researchers can write their own commands, packages are how they are distributed for other people to use
* `install.packages("dplyr")` you only need to install once per machine!
	* you may be asked to choose a site, but this doesn't matter much (RStudio mirror)
* `library("dplyr")` you'll need to load the package every time you reboot R
* may see red text output: these are probably not errors, just warnings!
* check to see if installation worked by typing `?select`

### What is dplyr?
* dplyr is a package that helps with common data manipulation tasks, especially data frames
* includes a number of functions, we'll learn a few of the most widely used

### Selecting columns and filtering rows
* `select(surveys, plot_id, species_id, weight)` for certain columns
* `filter(surveys, year == 1995)` for certain rows

### Pipes
* You may want to select rows and columns at the same time
* there are many ways to do this:
	* intermediate steps (create temp objects, may clutter up space)
	* nested functions (may be difficult to read)
	* pipes: take output from one function and send as input to another function
* pipes are a common programming model, and these data will help you become familiar with them
	* in R, pipes are `%>%` and are available through the package magrittr (installed with dplyr)
```
surveys %>%
filter(weight < 5) %>% 
select(species_id, sex, weight)
```
* survey data sent through filter, then that output is sent through select 
* can assign to object, view output
Challenge question

### Mutate
* create new columns based on existing columns
* create new column of weight in kg:
```
surveys %>%
	mutate(weight_kg = weight / 1000)
```
* pipe to head to see output
* insert filter to remove NAs:
```
surveys %>%
	filter(!is.na(weight)) %>% 
	mutate(weight_kg = weight / 1000) %>% 
	head
```
* **Challenge** Using pipes, subset the data to include rows before 1995. Retain columns year, sex, and weight.

### Split-apply-combine with summarize
* split data into groups, apply analysis to each group, combine results
* to count number of rows for each sex:
```
surveys %>% 
	group_by(sex) %>% 
	tally()
```
* collapse each group into single-row summary:
```
surveys %>%
	group_by(sex) %>%
	summarize(mean_weight = mean(weight, na.rm = TRUE)) #na.rm removes missing data
```
* group by multiple columns
```
surveys %>%
	group_by(sex, species_id) %>%
	summarize(mean_weight = mean(weight, na.rm = TRUE))
```
* some species weren't weighed, but we can filter out those too!
```
surveys %>%
	group_by(sex, species_id) %>%
	summarize(mean_weight = mean(weight, na.rm = TRUE)) %>% 
	filter(!is.nan(mean_weight))
```
* summarize multiple variables at the same time
```
surveys %>%
	group_by(sex, species_id) %>%
	summarize(mean_weight = mean(weight, na.rm = TRUE),
		min_weight = min(weight, na.rm = TRUE)) %>% 
	filter(!is.nan(mean_weight))
```
* Challenge: 
	* how many times was each plot_type surveyed?
	* Use group_by() and summarize() to find the mean, min, and max hindfoot length for each species.
	* What was the heaviest animal measured in each year? Return the columns year , genus , species , and weight
* dplyr cheatsheet: http://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf
	
## DATA VISUALIZATION WITH GGPLOT2
* **Objectives**: visualize mammals data, understand ggplot2 methods, build step-by-step complex plots with ggplot2

* load required packages
```
install.packages("ggplot2")
library("ggplot2")
```
* load data from figshare: `surveys_raw <- read.csv("https://ndownloader.figshare.com/files/2292172")`

### Data cleaning and preparing for plotting
* `summary(surveys_raw)`
* remove missing values for species ID:
```
surveys_complete <- surveys_raw %>% 
	filter(species_id != "")
```
* remove missing values for weight and hindfoot length
```
surveys_complete <- surveys_raw %>% 
	filter(species_id != "") %>% # remove missing species_id
	filter(!is.na(weight)) %>%  # remove missing weight
	filter(!is.na(hindfoot_length)) # remove missing hindfoot_length
```
* remove species with less than 10 counts
```
# count records per species 
species_counts <- surveys_complete %>% 
	group_by(species_id) %>%
	tally
head(species_counts)

# get names of frequent species
frequent_species <- species_counts %>% 
	filter(n >= 10) %>%
	select(species_id)
	
surveys_complete <- surveys_complete %>%
	filter(species_id %in% 
	frequent_species$species_id)
```
* make simple plot of hindfoot length as a function of weight
`plot(x = surveys_complete$weight, y = surveys_complete$hindfoot_length)`

### Plotting with ggplot2
* ggplot2 is a plotting package that helps create complex, publication quality plots with minimal effort
* add layers of complexity to show data in the desired manner
* bind plot to specific data frame:`ggplot(data = surveys_complete)`
* define aesthetics that map variables to axes: `ggplot(data = surveys_complete, aes(x = weight, y = hindfoot_length))`
* add geoms that represent data on the plot:`ggplot(data = surveys_complete, aes(x = weight, y = hindfoot_length)) + geom_point()`

### Modifying plots
* add transparency: `ggplot(data = surveys_complete, aes(x = weight, y = hindfoot_length)) + geom_point(alpha = 0.1)`
* add color: `ggplot(data = surveys_complete, aes(x = weight, y = hindfoot_length)) + geom_point(alpha = 0.1, color = "blue")`

### Boxplot
* visualize distribution of weight within each species: `ggplot(data = surveys_complete, aes(x = species_id, y = weight)) + geom_boxplot()`
* add points:
```
ggplot(data = surveys_complete, aes(x = species_id, y = weight)) +
	geom_jitter(alpha = 0.3, color = "tomato") +
	geom_boxplot(alpha = 0)
```
* **Challenges:**
	* Replace the box plot with a violin plot; see `geom_violin()`
	* Represent weight on the log10 scale; see `scale_y_log10()`
	* Create boxplot for hindfoot_length .
	
### Plotting time series data
* calculate number of counts per year for each species:
```
yearly_counts <- surveys_complete %>% 
	group_by(year, species_id) %>%
	tally
```
* timelapse data shows time on x axis and counts on y axis: `ggplot(data = yearly_counts, aes(x = year, y = n)) + geom_line()`
	* this shows all species together
* show species separately: `(data = yearly_counts, aes(x = year, y = n, group = species_id)) + geom_line()`
* add colors: `ggplot(data = yearly_counts, aes(x = year, y = n, group = species_id, color = species_id)) + geom_line()`

### Faceting
* split one plot into multiple plots based on a factor
* plot one time series for each species separately: `ggplot(data = yearly_counts, aes(x = year, y = n, color = species_id)) + geom_line() + facet_wrap(~species_id)`
* What if we wanted a separate line in each facet for male and female?
* **Challenges:**
	* filter the dataframe so that we only keep records with sex “F” or “M”s
	```
	sex_values = c("F", "M")
	surveys_complete <- surveys_complete %>%
		filter(sex %in% sex_values)
	```
	* group by year, species_id, sex
	```
	yearly_sex_counts <- surveys_complete %>%
		group_by(year, species_id, sex) %>%
		tally
	```
	* make faceted plit splitting further by sex (within single plot)
	```
	ggplot(data = yearly_sex_counts, aes(x = year, y = n, color = species_id, group = sex)) +
	geom_line() + facet_wrap(~ species_id)
	```
* change default background:
```
ggplot(data = yearly_sex_counts, aes(x = year, y = n, color = species_id, group = sex)) +
geom_line() + facet_wrap(~ species_id) + theme_bw()
```
* color by sex instead of species:
`ggplot(data = yearly_sex_counts, aes(x = year, y = n, color = sex, group = sex)) + geom_line() + facet_wrap(~ species_id) + theme_bw()`
* plot average weight of each species over the course of years
```
yearly_weight <- surveys_complete %>%
	group_by(year, species_id, sex) %>%
	summarise(avg_weight = mean(weight, na.rm = TRUE))
ggplot(data = yearly_weight, aes(x=year, y=avg_weight, color = species_id, group = species_id)) +
	geom_line() + theme_bw()
```
* lines are in steps because of multiple values by year
* make separate plots per sex since weight of males and females can differ a lot
```
ggplot(data = yearly_weight, aes(x=year, y=avg_weight, color = species_id, group = species_id)) +
	geom_line() + facet_wrap(~ sex) + theme_bw()
```
* other ways to improve:
	* https://www.rstudio.com/wp-content/uploads/2015/08/ggplot2-cheatsheet.pdf
* change axis names:
```
ggplot(data = yearly_sex_counts, aes(x = year, y = n, color = sex, group = sex)) + geom_line() +
	facet_wrap(~ species_id) +
	labs(title = 'Observed species in time',
		x = 'Year of observation',
		y = 'Number of species') + theme_bw()
```
* change size and font of text:
```
ggplot(data = yearly_sex_counts, aes(x = year, y = n, color = sex, group = sex)) + geom_line() +
	facet_wrap(~ species_id) +
	labs(title = 'Observed species in time',
		x = 'Year of observation',
		y = 'Number of species') 
	theme(text=element_text(size=16, family="Arial")) + theme_bw()
```
* change x-axis orientation
```
ggplot(data = yearly_sex_counts, aes(x = year, y = n, color = sex, group = sex)) + 	
	geom_line() +
	facet_wrap(~ species_id) +
	theme_bw() +
	theme(axis.text.x = element_text(colour="grey20", size=12, angle=90, hjust=.5, vj
ust=.5),
		axis.text.y = element_text(colour="grey20", size=12),
		text=element_text(size=16, family="Arial")) + 
	labs(title = 'Observed species in time',
          x = 'Year of observation',
          y = 'Number of species')
```
* save plot to file: `ggsave("observed_species_in_time.png", width=15, height=10)`
	* file format included in file name, width and height can also be specified

## END CLASS
* remove R history from dropbox (remind students to save first if they want it)
* give students links to complete R lessons (and intermediate lessons)
* sticky note feedback
