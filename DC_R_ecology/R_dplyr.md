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
