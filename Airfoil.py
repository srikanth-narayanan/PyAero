"""
This is module is to read an airfoil file and create an airfoil object. This 
object contains methods analyze the airfoil geometry.

Next version of this python module will contains methods to optimize the shape 
of the airfoil using different methods.

:date: Feb 13th 2015
"""

__author__ = "Srikanth Narayanan"
__version__ = "0.1.0"
__email__ = "srikanth.n.narayanan@gmail.com"

#import libraries
import os
import numpy as np

class Airfoil(object):
    '''
    This is airfoil class and is initialized by passing the path of an standard 
    airfoil foil. This class has methods
    '''
    def __init__(self, airfoil_file=None):
        '''
        Constructor to initialise Airfoil object.
        
        :param airfoil_file: a full absolute path of the airfoil file
        :type airfoil_file: a raw input path string
        '''
        if airfoil_file:
            self.airfoil_path=os.path.airfoil
            self.load_airfoil(self.airfoil_path)
        else:
            print "No Airfoil file provided"
    
    def load_airfoil(self,airfoil_path):
        '''
        Method to read the given airfoil file
        
        :param airfoil_path: a full absolute path of the airfoil file
        '''

        
if __name__ == "__main__":
    pass