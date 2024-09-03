import csv 
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.patches as patches
import matplotlib.transforms as transforms
import os
import imageio

# importing movie py libraries
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage

class simulationData():
    def __init__(self, filename):
        self.filename=filename
        self.timeStepData=[]
        
        #declaring here so I don't have to in getInfo method
        self.timeStepDict={}
    
    #Structure of dict
    """
    timestep={
        'label' = [x, y, angle, length, radius, damping]
    }
    """
  
    #method to sort data into dictionaries
    def getInfo(self):
        #getting the save time from params file
        paramPath='/'.join(self.filename.split('/')[:-1])
        with open(os.path.join(paramPath, f'params.txt'), 'r') as file:
            lines=file.read().strip().split('\n')

        # Extract headers and values
        headers=lines[0].split(', ')
        values=lines[1].split(', ')

        # Create a dictionary
        params_dict=dict(zip(headers, values))

        # Extract the save_time
        self.saveTime=float(params_dict['save_time'])
        
        # Open the file in read mode - getting data info
        with open(self.filename, 'r') as file:
            # Create a CSV reader object
            reader=csv.reader(file)

            self.numTimesteps=0
            
            # Iterate over each row in the reader
            for row in reader:
                # Check if the row is empty
                if len(row) == 0:
                    self.timeStepData.append(self.timeStepDict)
                    
                    #resetting dictionary for next timeframe
                    self.timeStepDict={}
                    self.numTimesteps+=1
                    
                else:
                    #put data for label into dictionary
                    self.timeStepDict[row[0]]=[float(row[1].strip()), #x
                                          float(row[2].strip()), #y
                                          float(row[3].strip()), #angle
                                          float(row[4].strip()), #length
                                          float(row[5].strip())] #radius
            
            #keeping track of the number of timesteps
            self.numTimesteps+=1
    
    #method to find parent/gparent/ggparent/gggparent
    def findClosestRelative(self, cellID, currentFrame): #current frame should always be at least 2, remember index will be different to currentFrame
        #checking if same cell in previous frame (not split yet)
        if cellID in self.timeStepData[currentFrame-1].keys(): 
            return (cellID, currentFrame-1)

        #checking if parent in previous frame (it has split in timestep)
        elif cellID[:-1] in self.timeStepData[currentFrame-1].keys():
            return (cellID[:-1], currentFrame-1)

        #checking if grandparent in previous frame
        elif cellID[:-2] in self.timeStepData[currentFrame-1].keys():
            return (cellID[:-2], currentFrame-1)
        
        #NOT SURE ABOUT THE FOLLOWING CASES - THEY ARE REQUIRED FOR CODE TO WORK WHEN RUNNING SIMULATION WHEN WE ONLY SAVE DATA EVERY SIMULATION HOUR
        #checking for great grandparent in double previous
        elif cellID[:-3] in self.timeStepData[currentFrame-2].keys():
            return (cellID[:-3], currentFrame-2)
        
        elif cellID[:-3] in self.timeStepData[currentFrame-1].keys():
            return (cellID[:-3], currentFrame-1)
        
        #gggparent
        elif cellID[:-4] in self.timeStepData[currentFrame-2].keys():
            return (cellID[:-4], currentFrame-2)
        
        elif cellID[:-4] in self.timeStepData[currentFrame-1].keys():
            return (cellID[:-4], currentFrame-1)
        
        elif cellID[:-4] in self.timeStepData[currentFrame-3].keys():
            return (cellID[:-4], currentFrame-3)
        
           
            
    #calculating Y RMS of single cell compared to its relative in each timestep    
    def singleCellRMS(self): #default range(2, len(self.timeStepData)) (for all timeframes)
        
        RMSStrain0Average=[] #for average value of RMS in each time frame
        RMSStrain1Average=[]
        
        for timestep in range(2, len(self.timeStepData)): #starting from frame 2, timestep=currentFrame
            
            RMSStrain0FrameValues=[] #for all values of RMS for current frame
            RMSStrain1FrameValues=[]
            
            #finding parents of each cell in a timestep
            for key in self.timeStepData[timestep].keys(): #key=cellID
                #getting currentInfo
                currentYInfo=self.timeStepData[timestep][key][1]
                
                
                #finding 'parent'
                parentInfo=self.findClosestRelative(key, timestep)
                parentYInfo=self.timeStepData[parentInfo[1]][parentInfo[0]][1]
                
                #finding 'grandparent'
                grandParentInfo=self.findClosestRelative(parentInfo[0], parentInfo[1])
                
                #handles for when we find parent in frame 1 and need to find grandparent that doesn't exist
                if grandParentInfo==None:
                    singleCellRMS=np.sqrt((currentYInfo-np.mean([parentYInfo]))**2)
                else:
                    grandParentYInfo=self.timeStepData[grandParentInfo[1]][grandParentInfo[0]][1]
                    singleCellRMS=np.sqrt((currentYInfo-np.mean([parentYInfo, grandParentYInfo]))**2)

                #keeping track of RMS of each strain
                if key[0]=='0':
                    RMSStrain0FrameValues.append(singleCellRMS)
                
                elif key[0]=='1':
                    RMSStrain1FrameValues.append(singleCellRMS)
                    
            #calculating the average
            RMSStrain0Average.append(np.mean(RMSStrain0FrameValues))
            RMSStrain1Average.append(np.mean(RMSStrain1FrameValues))
            
            
        return (RMSStrain0Average, RMSStrain1Average)
    
    #method to calculate whether we have coexistence or which strain won
    def calculateOutcome(self):
        strainOneSurvive=False #here one and zero refer to ID, not physical number 
        strainZeroSurvive=False
        
        #checking the last timestep for cell types
        for key in self.timeStepData[-1].keys():
            if key[0]=='1':
                strainOneSurvive=True
                break
            elif key[0]=='0':
                strainZeroSurvive=True
                break
         
        #categorising outcome
        if strainOneSurvive==True and strainZeroSurvive==True:
            self.result=2 #2 for coexist
        elif strainOneSurvive==True and strainZeroSurvive==False:
            self.result=1 
        elif strainOneSurvive==False and strainZeroSurvive==True:
            self.result=0
        
        return self.result
    
    #method to find the average ordering of each timestep
    def calculateOrdering(self):
        orderingAverage=[]
        
        #looking through each timestep
        for timestep in range(len(self.timeStepData)):
            angles=[]
            
            #looking at each cell in the timestep and getting its angle
            for key in self.timeStepData[timestep].keys():
                angles.append(self.timeStepData[timestep][key][2])
            
            #calculating ordering with equation 
            orderingAverage.append(np.average((3*np.cos(angles)**2-1)/(2)))
            
        return orderingAverage 
    
    #method to find ordering of individual strains in simulation
    def calculateStrainOrdering(self):
        str0orderingAverage=[]
        str1orderingAverage=[]
        
        for timestep in range(len(self.timeStepData)):
            str0angles=[]
            str1angles=[]
            
            for key in self.timeStepData[timestep].keys():
                #checking what type of bacteria we are looking at in loop
                if key[0]=='0':
                    str0angles.append(self.timeStepData[timestep][key][2])
                if key[0]=='1':
                    str1angles.append(self.timeStepData[timestep][key][2])
            
            str0orderingAverage.append(np.average((3*np.cos(str0angles)**2-1)/(2)))
            str1orderingAverage.append(np.average((3*np.cos(str1angles)**2-1)/(2)))
        
        return str0orderingAverage, str1orderingAverage
                
    
    #method to calculate the population of each strain at each timestep
    def calculatePopulations(self):
        pop0=[]
        pop1=[]
        
        #looking through timesteps
        for timestep in range(len(self.timeStepData)):
            num0=0; num1=0
            #counting number of each cell in the timestep
            for key in self.timeStepData[timestep].keys():
                if key[0]=='1':
                    num1+=1
                elif key[0]=='0':
                    num0+=1
            
            pop0.append(num0)
            pop1.append(num1)
        
        return (pop0, pop1)
    
    #method to get the point at which fixation occurs
    #this method currently assumes that timesteps are in intervals of hours
    def getFixationTime(self):
        self.calculateOutcome()
        pop0, pop1 = self.calculatePopulations()
        
        #no fixation time if cells coexist
        if self.result==2:
            self.fixationTime=None
        
        #if 0 wins, check when 1 dies
        if self.result==0:
            for i in range(len(pop1)):
                if pop1[i]==0:
                    self.fixationTime=i*(self.saveTime/60) #save time is in minutes in param file
                    break
        
        #if 1 wins, check when 0 dies
        if self.result==1:
            for i in range(len(pop0)):
                if pop0[i]==0:
                    self.fixationTime=i*(self.saveTime/60)
                    break
        
        return self.fixationTime
            
    #method to calculate fractional occupation of both strains at each timestep, quite a rough method
    def calculateFractionalOccupation(self):
        #calculating populations
        pop0, pop1 = self.calculatePopulations()
        
        #total population is equivalent to last timestep numbers added
        totalPopulation=pop0[-1]+pop1[-1]
        
        #convert to np array so we can divide by single number
        return np.array(pop0)/totalPopulation, np.array(pop1)/totalPopulation

  #alternative method to calculate fractional area, using a calculated constant and considering the simulation space  
  def calculateFractionalOccupationAlternative(self):
        totalArea=self.channelWidth*self.channelHeight
        
        pop0radii=[]; pop1radii=[]
        pop0length=[]; pop1length=[]
        
        #calculating populations
        pop0, pop1 = self.calculatePopulations()
        
        #getting average length and radii
        for key in self.timeStepData[-1]:

            if key[0]=='1':
                pop1radii.append(self.timeStepData[-1][key][4])
                pop1length.append(self.timeStepData[-1][key][3])
            elif key[0]=='0':
                pop0radii.append(self.timeStepData[-1][key][4])
                pop0length.append(self.timeStepData[-1][key][3])
        
        #getting average length and radii
        (pop0radiimu, pop0radiisigma) = scipy.stats.norm.fit(pop0radii)
        (pop0lengthmu, pop0lengthsigma) = scipy.stats.norm.fit(pop0length)
        
        (pop1radiimu, pop1radiisigma) = scipy.stats.norm.fit(pop1radii)
        (pop1lengthmu, pop1lengthsigma) = scipy.stats.norm.fit(pop1length)

        #calculating constant to multiply population size by
        pop0areaConst=(pop0lengthmu-2*pop0radiimu)*pop0radiimu+(np.pi*pop0radiimu**2)
        pop1areaConst=(pop1lengthmu-2*pop1radiimu)*pop1radiimu+(np.pi*pop1radiimu**2)
        
        #convert to np array so we can divide by single number
        return (np.array(pop0)*pop0areaConst)/totalArea, (np.array(pop1)*pop1areaConst)/totalArea
    

    #method to draw bacteria to be visualised
    def drawBacteria(self, ax, x, y, length, radius, angle, arrowBoolean, edgecolor='none', facecolor='none', arrow_length=0.5):
        # Calculate the end position of the arrow based on the angle
        arrow_dx=arrow_length*np.cos(np.radians(angle))
        arrow_dy=arrow_length*np.sin(np.radians(angle))

        #if bacteria is pill shaped
        if length>2*radius:
            #transforming length to be drawn correctly
            rectlength=length-(2*radius)

            rect=patches.Rectangle((x-(rectlength/2), y-(radius)), rectlength, 2*radius, angle=angle, rotation_point='center', facecolor=facecolor, edgecolor=edgecolor)
            ax.add_patch(rect) #fixed drawing issue by line above: x-rectlength rather than x-length


            left_cap=patches.Wedge((x-(rectlength/2*np.cos(np.radians(angle))), y-(rectlength/2*np.sin(np.radians(angle)))), radius, 90+angle, 270+angle, facecolor=facecolor, edgecolor=edgecolor)
            ax.add_patch(left_cap)

            right_cap=patches.Wedge((x+(rectlength/2*np.cos(np.radians(angle))), y+(rectlength/2*np.sin(np.radians(angle)))), radius, -90+angle, 90+angle, facecolor=facecolor, edgecolor=edgecolor)
            ax.add_patch(right_cap)
            #here radius=height, i.e i want the height of the circle wedge to be height of rectangle, hence height=0.5*rad of rectangle

        if length<2*radius:
            ellipse=patches.Ellipse((x, y), length, radius*2, angle=angle, edgecolor=edgecolor, facecolor=facecolor) #check this angle is right
            ax.add_patch(ellipse) #we want diameter here

        #if bacteria is circle shaped
        elif length==radius: #feel like i can use else statement here 
            circle=patches.Circle((x, y), radius / 2, edgecolor=edgecolor, facecolor=facecolor)
            ax.add_patch(circle)

        # Add the arrow to all shapes
        if arrowBoolean==True:
            ax.annotate('', xy=(x + arrow_dx, y + arrow_dy), xytext=(x, y), arrowprops=dict(facecolor='green', edgecolor='green', arrowstyle='->', lw=2))

    #method to draw frames of sim and save to folder
    def drawFrames(self, arrowBoolean):
        output_dir='frames'
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        
        frame_counter = 0
        
        fig, ax=plt.subplots(figsize=(10, 10))

        for timestep in range(len(self.timeStepData)):
            for key in self.timeStepData[timestep].keys():
            
                if key[0]=='1':
                    self.drawBacteria(ax, self.timeStepData[timestep][key][0], self.timeStepData[timestep][key][1], self.timeStepData[timestep][key][3], self.timeStepData[timestep][key][4], self.timeStepData[timestep][key][2], arrowBoolean, edgecolor='r', facecolor='r')
                if key[0]=='0':
                    self.drawBacteria(ax, self.timeStepData[timestep][key][0], self.timeStepData[timestep][key][1], self.timeStepData[timestep][key][3], self.timeStepData[timestep][key][4], self.timeStepData[timestep][key][2], arrowBoolean, edgecolor='b', facecolor='b')
                
            # Set the xlim and ylim for the current frame
            ax.set_xlim(0, 44)  # Set up size of space; this data should ideally be read from a parameter file
            ax.set_ylim(0, 12)  # as it could change
            ax.set_aspect('equal', adjustable='box')

            frame_filename=os.path.join(output_dir, f'frame_{frame_counter:04d}.png')
            plt.savefig(frame_filename)
            #print(f'Saved {frame_filename}')

            frame_counter+=1
            plt.clf()
            fig, ax=plt.subplots(figsize=(10, 10))   
        
        plt.close(fig)
        
    #method to make gif of frames of sim and save it to folder
    def createGif(self, fps):
        # Get a list of image file paths sorted by name
        image_files = sorted([os.path.join('frames', img) for img in os.listdir('frames') if img.endswith(".png")])

        # Read each image and store in a list
        images = []
        for filename in image_files:
            images.append(imageio.imread(filename))

        # Create the output directory if it doesn't exist
        if not os.path.exists('gif'):
            os.makedirs('gif')

        # Save the images as a GIF with the correct file path
        imageio.mimsave(os.path.join('gif', 'animation.gif'), images, fps=fps)
