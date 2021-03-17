## Online book https://r4ds.had.co.nz

## install and load tidyverse

# install.packages("tidyverse")
library(tidyverse)

## Take careful note of the conflicts message that’s printed when you load the tidyverse. 
## It tells you that dplyr overwrites some functions in base R. 
## If you want to use the base version of these functions after loading dplyr, 
## you’ll need to use their full names: stats::filter() and stats::lag().

################################################################################
## ggplot2
################################################################################

## First qplot

## look at the data:
## mpg (tibble!) details on variables
## https://ggplot2.tidyverse.org/reference/mpg.html 

qplot(x = displ, y = hwy,        # displ   car's engine size (liters),
      data = mpg)                # hwy    fuel efficiency on highway
                                 #        (low efficiency -> more fuel)
qplot(x = displ, y = hwy,
      data = mpg, geom = "line")
qplot(x = displ, y = hwy,
      color = drv, data = mpg)  # add colors for different classes,
                                #   f = front-wheel drive,
                                #   r = rear wheel drive,
                                #   4 = 4wd

qplot(x = displ,                # only x -> histogram
      data = mpg)         
qplot(x = displ,
      data = mpg, geom = "density")


## Now qqplot2

ggplot(data = mpg)              # nothing happens here / empty plot
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy))
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy, color = drv))

## Try different aesthetic (color, size, alpha, shape)
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy, alpha = drv))
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy, shape = drv))
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy, size = drv))

## Set global aesthetic
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy),
             color = "blue",
             shape = 23,
             fill = "red",
             stroke = 2) # check ?geom_point for available aesthetics
                         # check vignette("ggplot2-specs") for all of them

## Facets
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy)) +
  facet_wrap(~ drv, nrow = 1)

## Adding a smooth line (loess smoother)
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy)) +
  geom_smooth(mapping = aes(x = displ, y = hwy))
  
## works the same with aesthetics
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy, color = drv)) +
  geom_smooth(se = FALSE,
              mapping = aes(x = displ, y = hwy, linetype = drv, color = drv))

## more types of geom_" geom_bar, geom_boxplot, 
ggplot(data = mpg, mapping = aes(x = class, y = hwy)) +
  geom_boxplot()

## maps
nz <- map_data("nz") #check ?map_data
ggplot(data = nz, mapping = aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "white", color = "black")

## Changes axes, labels and add legends, themes
g <- ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy, color = drv))

## sizes
g + theme(text = element_text(size = 20))
g + theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 30),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 15))
## labels
g <- g + theme(text = element_text(size = 20)) +
        xlab("Engine size") +
        ylab("Efficiency") +
        ggtitle("Cars")

## legend (much more on position, background color is possible)
g + guides(color = guide_legend(title = "Drive")) +
    theme(legend.position = c(0.8,0.8))

## changing xlim, ylim
g + xlim(c(0, 8)) +
    ylim(c(0, 50)) +
    theme(legend.position = c(0.8, 0.8))

## change theme; see https://ggplot2.tidyverse.org/reference/ggtheme.html 
g + theme_minimal()
g + theme_light()
#for more see ?theme

################################################################################
## dplyr package
################################################################################

# install.packages("nycflights13")
library(nycflights13)
flights
## the data are saved in a "tibble"
## (int=integer, dbl=double, chr=character, dttm = data+tme,
##  lgl= T/F, fctr=factor)

## filter -> pick observations by their value
filter(flights, month == 1, day == 1)     # other options: >, <, !=, ==, <=, >=
filter(flights, month == 11 | month == 12)  # logical: | (or), & (and), ! (not)
# same as 
filter(flights, month %in% c(11, 12))

## reorder rows
arrange(flights, year, month, day)
arrange(flights, desc(arr_delay))

## select columns
select(flights, year, month, day)
## same as 
select(flights, year:day)
select(flights, -(year:day))

## add new variable
flights_sml <- select(flights, year:day, ends_with("delay"), distance, air_time)
mutate(flights_sml,
       gain  = arr_delay - dep_delay,
       speed = distance / air_time * 60)

## summarize 
summarize(flights, delay = mean(dep_delay, na.rm = TRUE))
by_month <- group_by(flights, month) 
summarize(by_month, delay = mean(dep_delay, na.rm = TRUE))

## summarize/group with pipe (pipe does not work with ggplot2, since this
##                            was written later)
## to make this plot from summary
by_dest <- group_by(flights, dest)
delay <- summarize(by_dest,
                   count = n(),
                   dist  = mean(distance, na.rm = TRUE),
                   delay = mean(arr_delay, na.rm = TRUE))

delay <- filter(delay, count > 20, dest != "HNL")

## the same with pipe %>%
delays <- flights %>%
  group_by(dest) %>% 
  summarize(count = n(),
            dist  = mean(distance, na.rm = TRUE),
            delay = mean(arr_delay, na.rm = TRUE) ) %>%
  filter(count > 20, dest != "HNL")

ggplot(data = delays, mapping = aes(x = dist, y = delay)) +
  geom_point(aes(size = count), alpha = 1/3) +
  geom_smooth(se = FALSE)


## tibbles: create, read, ...

## convert a data frame to tibble
## iris is a regular data frame
class(iris)
#convert it a tibble
iris_tl <- as_tibble(iris)
class(iris_tl)

## make a new tibble
tibble(x = 1:5, y = 1, z = x^2 + y)
tibble(':)'   = "smile",
       ' '    = "space",
       '2000' = "number",
       no     = 0)
## check all possible data formats (dates, strings, etc)

## by default 10 rows and not all columns (depends on the width of the screen)
## are printed; to customise:
flights %>% print(n = 5, width = Inf)
## to reset defaults: options(tibble.width = Inf) to print all columns or
## options(tibble.print_max = n, tibble.print_min = m) to print rows

## extracting variables
df <- tibble(x = runif(5), y = rnorm(5))

df$x
## same
df[["x"]]
## same 
df[[1]]
## same
df %>% .$x
df %>% .[[1]]

## read the data (package readr in tidyverse)
shhs <- read_tsv("shhs1.txt")
## compare to (different data types! integer/double)
as_tibble(read.table("shhs1.txt", header = T))
## check parse_* functions!

## tsv reads tab-delimeted, csv comma-delimted, etc
## options to omit first rows, not to treat first rows as col.names
read_tsv("shhs1.txt", skip = 1, col_names = FALSE)

## writing to a file
write_tsv(shhs, "shhs2.txt")

## there are more options to read the data from SPSS, STATA, SAS, Excel, etc
## check help

################################################################################
## tidy data (package tidyr in tidyverse)
################################################################################

## tidy data -> each variables has its own column,
##              each observation has its own row,
##              each value has its own cell

## The packages dplyr, ggplot2 are designed to work with tidy data.
## in many data sets: one variable is spread across multiple columns
## one observation is scattered across multiple rows;
## to solve: garter(), spread(), separate(), pull()

## spread() makes long tables shorter, gather() makes wide table narrower
table4a
## make a tidy tibble with country, year, cases using gather()
table4a %>% gather(`1999`, `2000`, key = "year", value = "cases")

## spread is the opposite to gather
table2
## each observation in country and year is spread across two rows
spread(table2, key = type, value = count)

## separate makes two variables from obe, wherever a separator characted apperas
table3
table3 %>% separate(rate,
                    into = c("cases", "population"),
                    sep = "/",
                    convert = TRUE)
## Note; without convert cases and population are characters!

## unite is the to opposite
table5
table5 %>% unite(new, century, year, sep = "")


## dealing with missing values
stocks <- tibble(year   = c(2015, 2015, 2015, 2015, 2016, 2016, 2016),
                 qtr    = c(1, 2, 3, 4, 2, 3, 4),
                 return = c(1.88, 0.59, 0.35, NA, 0.92, 0.17, 2.66))
stocks

## Note: a return in the 4th quarter of 2015 is explicitly missing (NA),
##       while return in the first quater of 2016 is implicitly missing
##       (not there at all)
## -> to make an implicit missing explicit:
stocks %>% complete(year, qtr)

## case study data set WHO: cases of Tuberculosis broken down by
## year, country, age, gender, diagnosis method
who

## country, iso2, iso3 are all the same specifying the country
## new_sp_* all seems to be values, not variables, so gather all these
## and count the cases:
who1 <- who %>% gather(new_sp_m014:newrel_f65,
                       key = "key",
                       value = "cases",
                       na.rm = TRUE) # remove missings, otherwise counts = NA
who1
## new_sp_m014: old or new case (all new),
##              type of TB (sp, rel, ep, sn),
##              sex (m,f) and age group (0-14, 15-24, etc)

## there are cases with newrel instead of new_rel: fix
who2 <- who1 %>% mutate(key = stringr::str_replace(key, "newrel", "new_rel"))
who2

## now separate variables
who3 <- who2 %>% separate(key, c("new", "type", "sexage"), sep = "_")
who3

## new, iso2, iso3 are redundant
who4 <- who3 %>% select(-new, -iso2, -iso3)
who4

## separate age and sex
who5 <- who4 %>% separate(sexage, c("sex", "age"), sep = 1)
who5

## same, all together:
who %>%
  gather(code, value, new_sp_m014:newrel_f65, na.rm = TRUE) %>%
  mutate(code = stringr::str_replace(code, "newrel", "new_rel")) %>%
  separate(code, c("new", "var", "sexage")) %>%
  select(-new, -iso2, -iso3) %>%
  separate(sexage, c("sex", "age"), sep = 1)


## merging datasets (inner, full, left, right joins)
## a tibble with three persons with age
base <- tibble(id  = 1:3,
               age = seq(55, 60, length = 3))
base

visits <- tibble(id      = c(rep(1:2, 2), 4),
                 visit   = c(rep(1:2, 2), 1),
                 outcome = rnorm(5))
visits
## Note: the third person did not show up, but
##       there is a 4th person, of whom we do not know the age.

## inner_join joins all the columns that have the same name
inner_join(base, visits)  # case sensitive: id != ID
## may specify
inner_join(base, visits, by = "id") 

## left join
left_join(base, visits)   # will have everything what is in base and
                          # if something is missing in visits - NA

## right joint
right_join(base, visits) #will have everything in visits, rest NA; same as left_join(visits, base) up to order of variables

## full join
full_join(base, visits)

## merging on several variables
hr_visits <- tibble(id    = rep(1:2, 2),
                    visit = rep(1:2, 2),
                    hr    = rpois(4, lambda = 100))

full_join(visits, hr_visits)
## Note: same as full_join(visits, hr_visits, by = c("id", "visits"))


## Suggested: work through the chapters strings, factors, data and times
##            in R for Data Science (Chapters 14-16)
