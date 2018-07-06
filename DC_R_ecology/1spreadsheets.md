## Data Carpentry Spreadsheets Instructor Crib Sheet

### Before class (instructor):
* student resources:
	* class website
	* etherpad

## INTRODUCTION

### Lesson objectives
* Good data entry practices - formatting data tables in spreadsheets
* How to avoid common formatting mistakes
* Approaches for handling dates in spreadsheets
* Basic quality control and data manipulation in spreadsheets
* Exporting data from spreadsheets

### Overview
* we can use spreadsheets for data entry, organization, subsetting/sorting, stats, plotting
* we will not cover stats, plotting, writing code (macros) in spreadsheet programs
* problems with spreadsheets: using for anything other than data entry causes problems with flexibility and reproducibility

## FORMATTING DATA TABLES IN SPREADSHEETS

* since we will be using code to analyze data, we need to think like a computer in how we record that data
* it's common to need to reorganize previously entered data

### Keeping track of your analyses
* create new tab in spreadsheet to record notes
* identify steps taken to organize and clean data

### Structuring data in spreadsheets
* variables in columns (what is being measured)
* observations in rows
* only one piece of info per cell
* leave raw data raw
* export cleaned data into text-based format like CSV (comma-separated values)

### Exercise: http://www.datacarpentry.org/spreadsheet-ecology-lesson/01-format-data/
* download data: https://ndownloader.figshare.com/files/2252083
* open and reorganize, keeping our guidelines in mind
* group discussion

## FORMATTING PROBLEMS

* Using multiple tables
* Using multiple tabs
* Not filling in zeros
* Using problematic null values
* Using formatting to convey information
* Using formatting to make the data sheet look pretty
* Placing comments or units in cells
* Entering more than one piece of information in a cell
* Using problematic field names
* Using special characters in data
* Inclusion of metadata in data table
* Date formatting

## DATES AS DATA

* can store dates in single column or multiple columns
* Excel sometimes does horrible things to dates

### Exercise: http://www.datacarpentry.org/spreadsheet-ecology-lesson/03-dates-as-data/
* dates tab, "Date collected" column
* Extract month, day year to new column using YEAR() MONTH() DAY()

### Preferred date format
* multiple options, ensure accuracy

### Data formats in spreadsheets
* dates stored as integers

### Exercise: export as CSV

### Advantages of alternative formatting
* YEAR, MONTH, DAY
* YEAR, DAY-OF-YEAR
* single string: YYYYMMDDhhmmss

## QUALITY CONTROL

### Quality assurance
* select cells/column to validate
* Data tab, Data Validation
* Allow box to select type of data
* include allowable values in Source

1. Select the plot_id column
2. On the Data tab select Data Validation
3. In the Allow box select Whole number
4. Set the minimum and maximum values to 1 and 24.

* try entering wrong value
* Input Message
* Error Alert, Style

### Quality control
* save raw data file or tab
* README file
* sorting: remember to expand!
* conditional formatting: Format, Conditional Formatting, 2-Color Scale
* these steps can also be done in R, OpenRefine, SQL

## EXPORTING DATA

* Excel format is .xls or .xlsx
	* proprietary
	* others can't open
	* different versions
	* journal, archive, granting agency requiremnets
	* same for OpenOffice
* to save as CSV
	* File, Save as
	* Format, "Comma Separated Values"
	* Check file name and location, Save
* can open csv files in Excel 
* Windows vs Unix line endings, save as Windows csv
* R can import xls, but with caveats
* don't HAVE to use commas as delimiters
