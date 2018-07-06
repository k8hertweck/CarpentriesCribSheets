## Data Carpentry Spreadsheets Instructor Crib Sheet

### Before class (instructor):
* student resources:
	* class website
	* etherpad

## INTRODUCTION

* OpenRefine, formerly Google Refine
* free, open source tool for dealing with messy data
* why OpenRefine
	* know what you did to data
	* actions reversible
	* doesn't modify original file
	* saves time for messy data
	* can repeat on other datasets
	* easily apply clustering algorithms
* using program: install, use Chrome, http://127.0.0.1:3333/ or http://localhost:3333
* website has introductory videos: http://openrefine.org

## WORKING WITH OPEN REFINE

* start program, launch
* can accept multiple file types, including web address
* download data: http://www.datacarpentry.org/OpenRefine-ecology-lesson/setup/
* data is modified version of portal data

1. Create Project, Get data from This Computer.
2. Choose Files, select Portal_rodents_19772002_scinameUUIDs.csv
3. Next>> under the browse button to upload the data
4. OpenRefine gives preview, Create Project>> (upper right)

### Faceting
* adding multiple filters
* look for errors in scientific name

1. scientificName column
2. Facet > Text facet
3. left panel shows all unique names
4. sort by name and count

### Clustering
* automated way to fix multiple errors

1. scientificName Text Fact, click Cluster
2. can change Method and Keying Function
3. key collision method and metaphone3 keying function, 3 cluster
4. Merge? then Merge Selected and Recluster

### Split
* make genus and species separate columns

1. scientificName column
2. arrow at top of scientificName, Edit Column > Split into several columns...
3. Separator box, replace the comma with a space
4. uncheck box to remove column

### Undo/Redo
* option on left side of screen

### Trim Leading and Trailing Whitespace
* words with spaces at beginning or end

1. column scientificName, choose Edit cells > Common transforms > Trim leading and trailing whitespace
2. Split step has disappeared
3. redo split

## Filtering and Sorting with OpenRefine

### Filtering
* extract a subset of data

1. Click the down arrow next to scientificName > Text filter. A scientificName facet will appear on the left margin.
2. Type in bai and press return. There are 48 matching rows of the original 35549 rows (and these rows are selected for the subsequent steps).
3. At the top, change the view to Show 50 rows. This way you will see all the matching rows.

#### Exercise: What scientific names (genus and species) are selected by this procedure? 
Do Facet > Text facet on the scientificName column after filtering. This will show that two names match your filter criteria. They are Baiomys taylori and Chaetodipus baileyi.
#### Exercise: How would you restrict this to one of the species selected?
To restrict to only one of these two species, you could make the search case sensitive or you could split the scientificName column into species and genus before filtering or you could include more letters in your filter.

### Excluding entries
* use drop-down menu > Facet > Text facet to create a new facet. Only the entries with names that agree with your Text filter will be included in this facet.
* faceting gives overview description, while filtering allows selection of that data to subset

### Sorting
* drop-down menu for the selected column shows Sort.... Select what you would like to sort by 

* can use on mulitple columns

## Examining Numbers in OpenRefine

### Numbers
* all values are text by default
* transform to other data types using Edit cells > Common transforms

### Numeric facet
* find possible errors in data entry (non-number values or blanks)

#### Exercise:
1. For a column you transformed to numbers, edit one or two cells, replacing the numbers with text (such as abc) or blank (no number or text).
2. Use the pulldown menu to apply a numeric facet to the column you edited. The facet will appear in the left panel.
3. Notice that there are several checkboxes in this facet: Numeric, Non-numeric, Blank, and Error. Below these are counts of the number of cells in each category. You should see checks for Non-numeric and Blank if you changed some values.
4. Experiment with checking or unchecking these boxes to select subsets of your data.

### Scatterplot facet
* use column pulldown menu to > Facet > Scatterplot facet

### Examine pair of columns
* click on square in Scatterplot Matrix

## Scripts from OpenRefine

## Saving and exporting a project

## Other resources
