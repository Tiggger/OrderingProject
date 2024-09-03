# OrderingProject
Code for analysing agent based simulations of bacteria growing in a microchannel - Summer Research Project at UTM

# How to use the simulationData class definition

You want to initialise the object, providing it with the directory where the simulation data is located. Until you run the getInfo() method, the simulation data will not be extracted and added to the object. The information is put into dictionaries, which makes data processing speedy. The format of the dictionary entries for each time frame is as follows: timestep={'label' = [x, y, angle, length, radius, damping]}, where each entry in the timestep dictionary is a cell in that timestep. There is a list which holds all of the dictionaries, which encompasses all of the data in that given simulation object. 

To process multiple simulations to calculate the average behaviour of these, you must iterate through the files in the directory where you save the simulation data. The code I wrote to do this looks like:

```python
import os
directory = 'directorytosimulation/data' #ensure that in this directory, the simulations are there and also the params.txt file
numSims = 301 

file_names = [os.path.join(directory, f'sim{i}.txt') for i in range(numSims)]
instances = [simulationData(file) for file in file_names]

#list to hold data of all simulations
orderings=[] #holds the ordering information of each simulation, the length of this list should be the same as numSims

for instance in instances:

    instance.getInfo()

    #calculate some value
    ordering=instance.calculateOrdering

#calculating the average of all 300 simulations
avgOrdering=np.average(orderings, axis=0) #ensure to add axis=0 to calculate the average of each of the 0th items, each of the 1st item in the list etc

# findClosestRelative

This method is never used explicitly, but is an essential function to allow for single cell y-RMS calculations. It covers all of the cases for finding relatives, more than would be required typically. But for the simulations I was saving the state of the simulation every simulation hour, so including all of these cases for finding a relative ensured that I was always able to find a relative to compare data with.

# singleCellRMS

A method to calculate the avarage y-RMS distance (of all cells at each time frame of the simulation, compared to the most recent locatable relative in a previous time frame (where ever we can find the first relative). The function returns the average y-RMS of each strain at each timestep (list of values, each value represents the average y-RMS) as a list.

# calculateOutcome

A method to look at the final time frame of the simulation and deduce the outcome of that simulation (strain 0 wins, strain 1 wins, coexistance). The function returns the result, where 0 means strain 0 won the simulation, 1 meanining strain 1 won and 2 meaning the outcome was coexistance. 

# calculateOrdering

Method looks at the angles of every cell at each timestep, and calculates the ordering of that timestep using the equation: np.average((3*np.cos(angles)**2-1)/(2)). The function returns the average ordering value at each timestep as a list.

# calculateStrainOrdering

This method does the same calculations as the calculateOrdering function, but it separates the data for each strain, so that you can look at the difference in average ordering of each strain in the simulation, rather than the average ordering of the strains together. It is good for capturing the difference in ordering of each cell type (circle vs pill for example). Returns two lists which have the average ordering of each strain at each timestep. 

# calculatePopulations

A method which counts the number of cells in each population at each timestep. Returns a list of the population at each timestep for each type of cell.

# getFixationTime

A method which calculates at which time fixation occurs. It checks the outcome of the simulation and also calculated the populations of each strain using the functions defined above. It checks if the outcome of the simulation is a fixation event first, and then if this is the case it looks for the first time frame in which the losing cell type disappears. It uses information provided in the params.txt file the simulations spit out to calculate the time in hours. The function returns a number, which is the fixation time.

# calculateFractionalOccupation

This method calculates the fractional occupation of each type of cell in the simulation at each timestep. It is a rough calculation as we are considering the fraction of the population at each timestep compared to the total population in the final timestep. It returns a list for each cell type with this fractional area at each timestep. 

# calculateFractionalOccupationAlternative

This is another way in which we can calculate the fractional occupation at each timestep which is a little more technical and will capture the behaviour better. It calculates the average radius and length of each cell type, and using these calculates the area of a cell which these average values. The fractional area is then calculated as the population of a given cell type at a time frame, multiplied by this average area constant, divided by the total area of the microchamber (which for the majority of my simulations was 44x12um). The method returns two lists, each of which has this more refined fractional area calculated at each timestep.

# drawBacteria

This method shouldn't ever be explicitly used, as it is used in the drawFrames() method. It just handles drawing the different shapes of bacteria (circular, elliptical, pill shaped). To draw the pills, I am drawing two half circles, connected by a rectangle. 

# drawFrames(arrowBoolean)

This method is used to draw each frame of the simulation data. arrowBoolean is a boolean input which controls whether to draw unit arrows pointing in the direction of the angle as defined in the data output of the simulation. The arrows all point in the same direction, I'm not sure this is correct but is how I have implemented this. The method will create a new directory within the current directory called 'frames', where it will save a png for each frame, and number them.

# createGif(fps)

This method will combine the frames that are created by calling the drawFrames(arrowBoolean) method. In order to create the gif, you need to first call drawFrames(arrowBoolean). The fps argument controls how many frames per second will be displayed by the gif. This method will create a directory within the directory called 'gif' and will save the gif here titled animation.gif.
