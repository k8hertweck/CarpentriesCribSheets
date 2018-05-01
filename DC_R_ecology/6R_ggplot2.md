## Data Carpentry R Ecology Instructor Crib Sheet

### Before class (instructor):
* set up RStudio
* set up Socrative questions 
	* teacher room: PAW5AYWM
  	* start quiz, teacher paced, disable student names
* put history in accessible place for students (e.g., Dropbox)
	
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
