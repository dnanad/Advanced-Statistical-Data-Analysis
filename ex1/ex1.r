library("tidyverse")
library("haven")
library("maptools")
library("raster")
library("rgdal")
library("gridExtra")

####(a)####
# The data childrenfinal.dta are given in the STATA format.Read the data into R using an appropriate function from the tidyverse.
childrendata <- read_dta("childrenfinal.dta") #method suitable for STATA data

# Next, remove all variables that start with \s\, \v" and \m", followed by a number (avoid listing all of them).
#remove all variables that start with "s", "v" and "m", followed by a number
childrendata_1<- childrendata %>% 
  dplyr::select(-matches("^[svm][0-9]"))


# Check all the remaining variables. Do all variables have reasonable variable type (character,factor, double, integer, etc)? Convert the variables to a suitable type, if necessary.

#convert labelled double into double variables
childrendata_2<- childrendata_1 %>%
  mutate_if(is.double, as.double)

### (b) ###
# Make a smaller tibble that contains variables hypage, ruralfacto, female, zstunt, zweight, zwast, adm2.
childrendata_3<- childrendata_2 %>%
  dplyr::select(c(hypage, ruralfacto, female, zstunt, zweight, zwast, adm2))


#Variable zstunt is the so-called Z-score for stunting and is defined as the height of a child standardised with the median and standard deviation of heights of children at the same age from a healthy population. Children with Z-score less than -2 are dfined to be stunted.
#Make a scatter plot of zstunt against hypage.Add a smooth line to the plot.
g1 <- ggplot(childrendata_3, aes(x = hypage, y = zstunt)) +
  geom_point(alpha=0.2) +
  geom_smooth(se = T, color = "blue") +
  geom_hline(yintercept=-2, linetype='twodash', color='red', size=1.3) +
  labs(x = "Age", y = "Z-score")

g1
ggsave('scatter.png', plot = g1)

# Comment on the results. 

#Make smooth plots of zstunt against hypage for females and males on one plot, add a suitable legend.
#Use different colors for males and females.
g2 <- ggplot(childrendata_3, aes(x = hypage, y = zstunt,  col = factor(female))) +
  geom_point(alpha = 0.3) +
  geom_smooth(se = F) +
  scale_colour_manual(name='gender',labels = c("male", "female"), values = c("orange", "blue")) +
  theme(legend.key = element_rect(fill = "white", colour = "black"))+
  geom_hline(yintercept=-2, linetype='twodash', color='red', size=1.3) +
  labs(x = "Age", y = "Z-score")

g2

#Plot zstunt against age for urban and rural children
g3 <- ggplot(childrendata_3, aes(x = hypage, y = zstunt,  colour = factor(ruralfacto))) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = F) +
  scale_colour_manual(name='area',labels = c("urban", "rural"), values = c("blue", "deeppink")) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) +
  geom_hline(yintercept=-2, linetype='twodash', color='red4', size=1.3) +
  labs(x = "Age", y = "Z-score")
g3

#Comment on the results.
#Experiment with different aesthetics, themes and font sizes for the plots, report your favourite(s).
g4<-grid.arrange(g2, g3,nrow=2, heights=c(12,12))
ggsave('scatter_1.png', plot = g4)




### (c) ### 
#Most of the following commands are taken from https://rpubs.com/spoonerf/countrymapggplot2



   
# Plot the map of Kenya with all counties listed in adm2. Colour the county areas according to the mean of zstunt in the corresponding county.#Kenya shapefile data
Kenya<-getData("GADM", country="KE", level=1) 
plot(Kenya)#basic map

# setting an appropriate projection
Kenya_UTM<-spTransform(Kenya, CRS("+init=epsg:32537"))  

#the names of the regions
Kenya_UTM@data$NAME_1
#childrendata_3[7]
#sort the data in alphabetic order with respect to the column
childrendata_3<- childrendata_3[order(childrendata_3$adm2),]
Kenya_UTM@data<- Kenya_UTM@data[order(Kenya_UTM@data$NAME_1),]

#summarising childrendata_3 data according to the counties
childrendata_4 <- childrendata_3 %>% group_by(adm2) %>%
  summarise(Mzstunt = mean(zstunt,na.rm=TRUE), n = n())

#Note that the one county (Isiolo) is missing in the data.
#adding the missing county Isiolo
childrendata_4[nrow(childrendata_4) + 1,] <- NaN
childrendata_4$adm2[47] <- "ISIOLO"

#again reordering 
childrendata_4<- childrendata_4[order(childrendata_4$adm2),]

#dataframe for ggplot
Kenya_UTM@data$id <- rownames(Kenya_UTM@data)
Kenya_UTM@data <- mutate(Kenya_UTM@data, Mzstunt= childrendata_4$Mzstunt)
Kenya_df <- fortify(Kenya_UTM)
Kenya_df <- full_join(Kenya_df,Kenya_UTM@data, by="id")

#Make suitable legend and add county names (or corresponding labels) to the map. Comment on the results. In which counties are the children stunted?
# lets try to draw the map with colors according to zscore

ggplot() + 
  geom_polygon(data = Kenya_df, aes(x = long, y = lat, group = group, fill =
                                      Mzstunt), color = "black", size = 0.25) +
  theme(aspect.ratio = 1)+
  scale_fill_gradient(name='mean Z-stunt', high='white', low='darkred', na.value='grey50')+
  theme_void()+
  theme(aspect.ratio = 1)
  

#now we add names of the counties
#In order to add names to map, we need another dataframe with all the conunties' centroids
# "coordinates" extracts centroids of the polygons, in the order listed at Kenya1_UTM@data
centroids_df <- as.data.frame(coordinates(Kenya_UTM))
names(centroids_df) <- c("long", "lat")
childrendata_4<- childrendata_4[order(childrendata_4$adm2),]
centroids_df$NAME_1 <- Kenya_UTM@data$NAME_1
centroids_df$Mzstunt <- childrendata_4$Mzstunt

#Generating the map
ggplot(data = Kenya_df, aes(x = long, y = lat, group = group, fill = Mzstunt)) + 
  geom_polygon(color = "black", size = 0.25) +
  geom_text(data = centroids_df, aes(x = long, y = lat, label = NAME_1, group = NULL), size = 3) +
  scale_fill_distiller(name='mean Z-stunt', palette = "Spectral") +
  #scale_fill_gradient(name='mean Z-stunt', high='white', low='darkred', na.value='grey60')+
  theme_void() +
  theme(aspect.ratio = 1)

ggsave('Kenya.png')


##(d)write the tibble from (b) into a text file
write.table(childrendata_3,"childrendata_ex9.txt")
write_csv2(childrendata_3, 'children_ex9.csv')



