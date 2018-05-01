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
