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
* default data structure for tabular data, used for statistics and plotting
* import using `read.csv()` or `read.table()`; many default options
* import to know how R is interpreting your data!
* recall with `surveys`
* limited output: `head(surveys)`
* data.frame: table where columns are vectors of same length
* Object class:
	* `class(surveys)` - information about the class of the object
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
* **Challenge question** Based on the output of `str(surveys)`:
	* What is the class of the object `surveys`?
	* How many rows and how many columns are in this object?
	* How many species have been recorded during these surveys?

### Indexing and subsetting data frames
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
* **Challenge question** Create a data.frame containing only the observations from row 200 of the surveys dataset.

### Factors
* factors represent categorical data, can be ordered or unordered
* contain pre-defined set of values known as levels (assigned alphabetically)
* `sex <- factor(c("male", "female", "female", "male"))`
* `levels(sex)` and `nlevels(sex)`
* may need to specify the levels
* `sex <- factor(sex, levels = c("male", "female"))`
* converting factors
* renaming factors
* using `stringsAsFactors=FALSE`

### Formatting dates
* requires lubridate package
