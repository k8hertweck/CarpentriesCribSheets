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