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
library("tidyverse")
```
* load data: `surveys_complete <- read_csv("data_output/surveys_complete.csv")`

### Plotting with ggplot2
* make simple plot of hindfoot length as a function of weight
`plot(x = surveys_complete$weight, y = surveys_complete$hindfoot_length)`
* ggplot2 is a plotting package that helps create complex, publication quality plots with minimal effort
* add layers of complexity to show data in the desired manner
* bind plot to specific data frame:`ggplot(data = surveys_complete)`
* define aesthetics that map variables to axes: `ggplot(data = surveys_complete, aes(x = weight, y = hindfoot_length))`
* add geoms that represent data on the plot:`ggplot(data = surveys_complete, aes(x = weight, y = hindfoot_length)) + geom_point()`

### Modifying plots
* add transparency: `ggplot(data = surveys_complete, aes(x = weight, y = hindfoot_length)) + geom_point(alpha = 0.1)`
* assign plot to a variable: `surveys_plot <- ggplot(data = surveys_complete, aes(x = weight, y = hindfoot_length))`
* draw plot: placement of + matters!
```
surveys_plot + 
	geom_point()
```

### Building plots iteratively
* define datasets, lay out axes, choose geom:
```
ggplot(data = surveys_complete, aes(x = weight, y = hindfoot_length)) +
    geom_point()
```
* add transparency:
```
ggplot(data = surveys_complete, aes(x = weight, y = hindfoot_length)) +
    geom_point(alpha = 0.1)
```
* add color: 
```
ggplot(data = surveys_complete, aes(x = weight, y = hindfoot_length)) + 
	geom_point(alpha = 0.1, color = "blue")`
```
* color by species:
```
ggplot(data = surveys_complete, aes(x = weight, y = hindfoot_length)) +
    geom_point(alpha = 0.1, aes(color = species_id))
```
* Challenge: Use what you just learned to create a scatter plot of weight over species_id with the plot types showing in different colors. Is this a good way to show this type of data?

### Boxplot
* visualize distribution of weight within each species: `ggplot(data = surveys_complete, aes(x = species_id, y = weight)) + geom_boxplot()`
* add points:
```
ggplot(data = surveys_complete, aes(x = species_id, y = weight)) +
	geom_boxplot(alpha = 0) +
	geom_jitter(alpha = 0.3, color = "tomato")
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
* split the line in each plot by sex:
```
yearly_sex_counts <- surveys_complete %>%
                      group_by(year, species_id, sex) %>%
                      tally()
```
* make faceted plot:
```
ggplot(data = yearly_sex_counts, aes(x = year, y = n, color = sex)) +
     geom_line() +
     facet_wrap(~ species_id)
```
* change background:
```
ggplot(data = yearly_sex_counts, aes(x = year, y = n, color = sex)) +
     geom_line() +
     facet_wrap(~ species_id) +
     theme_bw() +
     theme(panel.grid = element_blank())
```

### Customization
* compare how weights of males and females have changed over time:
```
# One column, facet by rows
yearly_sex_weight <- surveys_complete %>%
    group_by(year, sex, species_id) %>%
    summarize(avg_weight = mean(weight))
ggplot(data = yearly_sex_weight, aes(x = year, y = avg_weight, color = species_id)) +
    geom_line() +
    facet_grid(sex ~ .)
```
* one row, facet by column:
```
# One row, facet by column
ggplot(data = yearly_sex_weight, aes(x = year, y = avg_weight, color = species_id)) +
    geom_line() +
    facet_grid(. ~ sex)
```
* change axes name:
```
ggplot(data = yearly_sex_counts, aes(x = year, y = n, color = sex)) +
    geom_line() +
    facet_wrap(~ species_id) +
    labs(title = "Observed species in time",
         x = "Year of observation",
         y = "Number of species") +
    theme_bw()
```
* increase font size:
```
ggplot(data = yearly_sex_counts, aes(x = year, y = n, color = sex)) +
    geom_line() +
    facet_wrap(~ species_id) +
    labs(title = "Observed species in time",
        x = "Year of observation",
        y = "Number of species") +
    theme_bw() +
    theme(text=element_text(size = 16))
```
* change axes labels to 90 degree angles
```
ggplot(data = yearly_sex_counts, aes(x = year, y = n, color = sex)) +
    geom_line() +
    facet_wrap(~ species_id) +
    labs(title = "Observed species in time",
        x = "Year of observation",
        y = "Number of species") +
    theme_bw() +
    theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
                        axis.text.y = element_text(colour = "grey20", size = 12),
          text = element_text(size = 16))
```
* can also change fonts and save themes
* save plot: `ggsave("fig_output/name_of_file.png", my_plot, width = 15, height = 10)`
* other ways to improve:
	* https://www.rstudio.com/wp-content/uploads/2015/08/ggplot2-cheatsheet.pdf

## END CLASS
* remove R history from dropbox (remind students to save first if they want it)
* give students links to complete R lessons (and intermediate lessons)
* sticky note feedback
