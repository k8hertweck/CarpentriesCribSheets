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
* tidyverse is a collection of packages that are very popular for data manipulation
* `install.packages("tidyverse")` you only need to install once per machine!
	* you may be asked to choose a site, but this doesn't matter much (RStudio mirror)
* `library("tidyverse")` you'll need to load the package every time you reboot R
* may see red text output: these are probably not errors, just warnings!
* check to see if installation worked by typing `?select`
* load data: `surveys <- read_csv("data/portal_data_joined.csv")`

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
	filter(!is.na(weight)) %>%
	group_by(sex, species_id) %>%
	summarize(mean_weight = mean(weight))
```
* add `print(n=15)` to end of last pipe to only print first 15 rows
* summarize multiple variables at the same time
```
surveys %>%
	filter(!is.na(weight) %>%
	group_by(sex, species_id) %>%
	summarize(mean_weight = mean(weight),
		min_weight = min(weight))
```
* add `arrange(min_weight)` to end of last pipe to order results
* count: `surveys %>%
			count(sex, sort=TRUE)`
* Challenge: 
	* how many times was each plot_type surveyed?
	* Use group_by() and summarize() to find the mean, min, and max hindfoot length for each species.
	* What was the heaviest animal measured in each year? Return the columns year , genus , species , and weight
* dplyr cheatsheet: http://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf

### Reshaping with gather and spread
* tie in with spreadsheet lesson
* use of tidyr package (installed with tidyverse)
* spreading will reshape the data (move rows to columns)
* create new object with mean weight:
```surveys_gw <- surveys %>%
 	 filter(!is.na(weight)) %>%
 	 group_by(genus, plot_id) %>%
 	 summarize(mean_weight = mean(weight))
```
* spread data:
```
surveys_spread <- surveys_gw %>%
  spread(key = genus, value = mean_weight)
str(surveys_spread)
```
* fill in missing values:
```
surveys_gw %>%
  spread(genus, mean_weight, fill = 0) %>%
  head()
```
* gather: want to treat columns as rows
```
surveys_gather <- surveys_spread %>%
  gather(key = genus, value = mean_weight, -plot_id)
```
* specify which columns to include:
```
surveys_spread %>%
  gather(key = genus, value = mean_weight, Baiomys:Spermophilus) %>%
  head()
```
* Challenge: Spread the surveys data frame with year as columns, plot id as rows, and the number of genera per plot as the values. You will need to summarize before reshaping, and use the function `n_distinct()` to get the number of unique genera within a particular chunk of data. 
```
rich_time <- surveys %>%
  group_by(plot_id, year) %>%
  summarize(n_genera = n_distinct(genus)) %>%
  spread(year, n_genera)
head(rich_time)
```

### Exporting data
* this also create the set of complete cases for our data visualization
```
surveys_complete <- surveys %>%
  filter(!is.na(weight),           # remove missing weight
         !is.na(hindfoot_length),  # remove missing hindfoot_length
         !is.na(sex))                # remove missing sex
```
* extract only most common species:
```
## Extract the most common species_id
species_counts <- surveys_complete %>%
    count(species_id) %>% 
    filter(n >= 50)

## Only keep the most common species
surveys_complete <- surveys_complete %>%
  filter(species_id %in% species_counts$species_id)
```
* save data to file: `write_csv(surveys_complete, path = "data_output/surveys_complete.csv")`