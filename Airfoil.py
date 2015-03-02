#!/usr/bin/env python
'''
This module contains class for creating Airfoil objects that can be used as
part of the aerodynamic framework. The objects of the Airfoil class support in
parsing a regular airfoil unit chord co-ordinate files.

:date: Feb 1, 2015
'''
__author__ = "Srikanth Narayanan"
__version__ = "0.2.1"
__email__ = "srikanth.n.narayanan@gmail.com"

import numpy as np
import os
from scipy.interpolate import interp1d as _i1d
from scipy.interpolate import InterpolatedUnivariateSpline as _IUS
from scipy.interpolate import UnivariateSpline as _US
import matplotlib.pyplot as plt
from matplotlib import rcParams as rcp
from scipy import stats
import math


class Airfoil(object):
    '''
    This class generates an airfoil object by reading a standard airfoil file.
    A standard airfoil file follows the trailing edge co-ordinate of the
    suction side, suction side leading edge, pressure side leading edge and
    pressure side trailing edge.

    The object also checks for airfoils in reversed direction and re-arrange
    the array to a default setting.
    Other methods included to calculate the following properties
    1. Trailing edge angle
    2. Maximum thickness
    3. Pressure side thickness
    4. Suction side thickness
    5. Pressure side thickness X-Co-ordinate.
    6. Suction side thcikenss Y-Co-ordinate.
    7. Leading edge radius.
    8. Airfoil area.
    9. Pressure side Inflection point.
    10. Airfoil Camber line.
    11. Trailing edge thickness.

    This class also support plotting the airfoil and save them to given a user
    directory as a png image file.
    '''

    def __init__(self, airfoil_file=None, logger=None, scons_pass=None):
        '''
        Constructor to initialise the input file.

        :param logger: Object of an airfoil logger.
        :param airfoil_file: An absolute path of the airfoil file.
        :param scons_pass: A flag to skip the logger initiation if called from
        sconstruct framework.
        '''
        if scons_pass is None:
            # setting up logger for GUI
            if logger:
                self.foil_io_logger = logger
            else:
                print "AIRFOIL.IO : No logger object provided"
            if airfoil_file:
                self.airfoil_file = airfoil_file
                self.load_airfoilfile(self.airfoil_file)
            else:
                print "AIRFOIL.IO : No input file given"

    def load_airfoilfile(self, airfoil_file):
        '''
        This method is to parse the Airfoil input file into a numpy array.
        Methods check for reverserd airfoil. It calculates trailing edge
        thickness.

        :param airfoil_file: The airfoil absolute file path to be read.
        '''
        self.airfoil_file = airfoil_file
        try:
            with open(airfoil_file, 'r') as fid:
                # read headers
                fid.readline()
                # read until end
                airfoil_table = []
                line = fid.readline()
                while line:
                    try:
                        if line != "\n":
                            data_line = [float(val) for val in line.split()]
                            airfoil_table.append(data_line)
                    except ValueError:
                        self._myprint("Error reading Airfoil profile file. \
Conversion error to float")
                        break
                    # continue reading
                    line = fid.readline()
                # convert table to numpy array
                self.airfoil_table = np.array(airfoil_table)
                # check and flip if airfoil table in reverse order
                check_status = self._check_foil_reverse()
                if check_status:
                    self._myprint("Airfoil in reverse order. Will be reversed \
for correct calcualtion of parameters")
                    self.airfoil_table = np.flipud(self.airfoil_table)
                self._myprint("Read airfoil file {0} successfully"
                              .format(airfoil_file))
        except:
            self._myprint("Error reading Airfoil input file")

    def _check_foil_reverse(self):
        '''
        This is a helper method to check if the airfoil is in reverse order.
        Forward order = x,y starting from suction side trailing edge and ending
        at pressure side trailing edge.

        A simple algorithm to check if the airfoil is in reverse order.
        Check the number of Y - cordinates in negative y-axis.

        Calculate the sum of number of negative point in the first half and
        second half and if the sum of first half is greater than second the
        airfoil is reversed.
        '''
        temp_y = self.airfoil_table[:, 1]
        con_mat = temp_y < 0
        n_points = int(len(temp_y)/2)
        sum_1sthalf = np.sum(con_mat[0:n_points])
        sum_2ndhalf = np.sum(con_mat[n_points:])

        # check condition
        if sum_1sthalf > sum_2ndhalf:
            return True
        else:
            return False

    def write_airfoil(self, fname, airfoil_cord, foil_header=None):
        '''
        This method writes the airfoil co-ordinates in to a given file name.

        :param fname: Absolute path and of the airfoil file to be written
        :param airfoil_cord: A numpy array of airfoil co-ordinates
        :param foil_header: Airfoil name as string variable.
        '''
        try:
            with open(fname, 'r') as fid:
                if foil_header is not None:
                    fid.write(foil_header)
                else:
                    fid.write('\n')
                np.savetxt(fid, airfoil_cord, fmt='%0.4f', delimiter='    ')
        except:
            self.foil_io_logger.error("Unable to write given airfoil \
co-ordinates")

    def find_trail_angle(self):
        '''
        This method will find the trailing edge angle. This method find the
        tangent equation to the trailing edge point and find the intercept
        angle between them.
        '''
        # study suggest use 5% of the total points
        npoints = len(self.airfoil_table)*5/100

        # find the slope of line1 close to trailing edge
        line1_m1, i1, r1, p1, std1 = stats.linregress(self.airfoil_table
                                                      [0:npoints, :])

        # find the slope of line2 formed by last three points close to
        # trailing edge
        line2_m2, i2, r2, p2, std2 = stats.linregress(self.airfoil_table
                                                      [-npoints:, :])
        self.trail_angle = np.abs(np.rad2deg(np.arctan((line1_m1 - line2_m2) /
                                                       (1 + (line1_m1 *
                                                        line2_m2)))))

    def thick_find(self):
        '''
        This method will find the thickness distribution of the airfoil.
        '''
        # find the index that seperates pressure side and suction side
        table_length = len(self.airfoil_table)
        # check if there even or odd number of points
        if table_length % 2 == 0:
            suc_side_idx = int(len(self.airfoil_table)/2)
            count_status = "Even"
        else:
            suc_side_idx = int(len(self.airfoil_table)/2)
            count_status = "Odd"
        # set pressure side start index
        if count_status == "Even":
            pre_side_idx = suc_side_idx
        if count_status == "Odd":
            pre_side_idx = suc_side_idx + 1

        # seperate data table
        self.X = self.airfoil_table[:, 0]
        self.Y = self.airfoil_table[:, 1]
        # find pressure and suction side x and y
        self.suc_side_x = self.X[:suc_side_idx]
        self.suc_side_y = self.Y[:suc_side_idx]
        self.pre_side_x = np.flipud(self.X[pre_side_idx:])
        self.pre_side_y = np.flipud(self.Y[pre_side_idx:])

        # find total thickness between pressure and suciton side, max thickness
        # and postion
        suc_spline = _IUS(np.flipud(self.suc_side_x),
                          np.flipud(self.suc_side_y))
        pre_spline = _IUS(np.flipud(self.pre_side_x),
                          np.flipud(self.pre_side_y))
        self.even_x = np.arange(0, 1, 0.005)
        suc_y = suc_spline(self.even_x)
        pre_y = pre_spline(self.even_x)
        self.thick = np.power((np.power((self.even_x - self.even_x), 2) +
                               np.power((suc_y - pre_y), 2)), 0.5)
        self.max_thick = np.nanmax(self.thick)
        self.max_thick_idx = np.nanargmax(self.thick)
        self.max_thick_x = self.suc_side_x[self.max_thick_idx]

        # create a unit chord and find pressure side and suction thickness
        self.chord_y = np.zeros(np.shape(self.suc_side_x))
        self.suc_side_thick = np.power((np.power((self.suc_side_x -
                                                  self.suc_side_x), 2) +
                                        np.power((self.suc_side_y -
                                                  self.chord_y), 2)), 0.5)
        self.pre_side_thick = np.power((np.power((self.pre_side_x -
                                                  self.pre_side_x), 2) +
                                        np.power((self.pre_side_y -
                                                  self.chord_y), 2)), 0.5)
        self.max_thick_suc = np.amax(self.suc_side_thick)
        self.max_thick_pre = np.amax(self.pre_side_thick)
        self.max_thick_suc_x = self.suc_side_x[self.suc_side_thick.argmax()]
        self.max_thick_pre_x = self.pre_side_x[self.pre_side_thick.argmax()]

    def find_pre_side_infl(self):
        '''
        This method finds the pressure side inflection point for the airfoil.
        '''
        # find pressure side point of inflection
        if hasattr(self, 'pre_side_x'):
            try:
                self.pre_infl_point = self._inflection_point(
                                      np.flipud(self.pre_side_x),
                                      np.flipud(self.pre_side_y))
            except Exception, err:
                self._myprint(str(err))
        else:
            self.thick()
            try:
                self.pre_infl_point = self._inflection_point(
                                      np.flipud(self.pre_side_x),
                                      np.flipud(self.pre_side_y))
            except Exception, err:
                self._myprint(str(err))

    def find_x_at_thick(self, percent_chord):
        '''
        This method find thickness at given percentage of chord.

        :param percent_chord: percentage chord as user input
        '''
        if hasattr(self, 'thick'):
            thick_ipl = _i1d(self.even_x, self.thick)
            return float(thick_ipl(percent_chord/100))
        else:
            return np.nan

    def find_LE_radius(self):
        '''
        This method finds the radius of the leading edge of the airfoil.
        Find the curvature of each point and inverse it to get radius of
        curvature.
        '''
        count = len(self.airfoil_table)/2
        try:
            xx = self.X[count-5:count+5]
            yy = self.Y[count-5:count+5]
            dx = np.gradient(xx)
            dy = np.gradient(yy)
            ddx = np.gradient(dx)
            ddy = np.gradient(dy)
            num = dx*ddy - ddx*dy
            den = (dx*dx)+(dy*dy)
            den = np.power(den, 1.5)
            self.curv_list = num/den
            self.le_radius = 1/np.amax(self.curv_list)
            self.le_point_idx = self.curv_list.argmax()
        except Exception, err:
            self._myprint(str(err))
            self.le_radius = np.nan

    def find_airfoil_area(self):
        '''
        This method find the area enclosed by airfoil.
        Integrate the area enclosed by suction side curve and pressure side
        curve seperately and add them together to get the total area.
        '''
        try:
            suc_side_spline = _IUS(np.flipud(self.suc_side_x),
                                   np.flipud(self.suc_side_y))
            pre_side_spline = _IUS(np.flipud(self.pre_side_x),
                                   np.flipud(self.pre_side_y))
            area1 = suc_side_spline.integral(0, 1)
            area2 = pre_side_spline.integral(0, 1)
            if not(math.isnan(area1) or math.isnan(area2)):
                self.airfoil_area = np.abs(area1-area2)
            else:
                # redundancy if interpolated spline fails use univariate spline
                suc_side_spline = _US(np.flipud(self.suc_side_x),
                                      np.flipud(self.suc_side_y))
                pre_side_spline = _US(np.flipud(self.pre_side_x),
                                      np.flipud(self.pre_side_y))
                area1 = suc_side_spline.integral(0, 1)
                area2 = pre_side_spline.integral(0, 1)
                if area1 != np.nan and area1 != np.nan:
                    self.airfoil_area = np.abs(area1 - area2)
                else:
                    self.airfoil_area = np.nan
        except Exception, err:
            self._myprint(str(err))
            self.airfoil_area = np.nan

    def find_stiffness_metrics(self, prut_chord_per):
        '''
        This method find the stifness parameter when the user defines the
        carbon slab postion as percentage of the chord length.

        :param prut_chord_per: Carbon slab width as a percentage chord
        '''
        try:
            set_inter = _i1d(self.even_x, self.thick)
            min_x = self.even_x[self.max_thick_idx] - ((prut_chord_per/100)/2)
            max_x = self.even_x[self.max_thick_idx] + ((prut_chord_per/100)/2)
            rnge_x = np.linspace(min_x, max_x, 10)
            ip_thick = set_inter(rnge_x)
            numera = (np.sum(np.power(ip_thick, 2)))/len(ip_thick)
            self.flap_stiff_metric = (np.power(numera, 0.5))
        except Exception, err:
            self._myprint(str(err))

    def _inflection_point(self, x, y):
        '''
        This method finds the point of inflection for the given x and y
        co-ordinate.

        :param x: the x axis co-ordinate for the given curve
        :param y: the y axis co-ordinate for the given curve
        '''
        try:
            try:
                inter_sp = _IUS(x, y)
                roots_ips = inter_sp.roots()
            except:
                # defensive when Interpolate Univariate Spline fails use
                # Univariate Spline
                inter_sp = _US(x, y)
                roots_ips = inter_sp.roots()
            tole = 0.1
            if len(roots_ips) > 0:
                if len(roots_ips) > 1:
                    new_roots = []
                    for val in roots_ips:
                        if val > tole:
                            new_roots.append(val)
                    return new_roots[0]
                else:
                    return roots_ips[0]
            else:
                # method to find inflection point if the spline does not touch
                # x-axis and doesn't have a root.
                samp_x = np.linspace(0, 1, 500)
                samp_y = inter_sp(samp_x)
                roots = samp_x[samp_y.argmax()]
                return roots
        except Exception, err:
            self._myprint(str(err))
            return np.nan

    def find_camber(self):
        '''
        This method will find the camber points for a given airfoil.
        '''
        if hasattr(self, "suc_side_x"):
            self._help_camber()
        else:
            self.thick_find()
            self._help_camber()

    def _help_camber(self):
        '''
        Camber helper method
        '''
        self.camber_x = (self.suc_side_x + self.pre_side_x)/2
        self.camber_y = (self.suc_side_y + self.pre_side_y)/2
        self.max_camber = np.amax(self.camber_y)
        self.max_camber_x = self.camber_x[self.camber_y.argmax()]
        self.camber_angle_TE = self._find_angle(self.camber_x[0:4],
                                                self.camber_y[0:4],
                                                self.camber_x[0:4],
                                                np.array([0, 0, 0, 0]))
        self.camber_angle_LE = self._find_angle(self.camber_x[-4:, ],
                                                self.camber_y[-4:, ],
                                                self.camber_x[-4:, ],
                                                np.array([0, 0, 0, 0]))

    def _find_angle(self, x1, y1, x2, y2):
        '''
        This is a helper method to find the angle between two sets of lines.

        :param x1: X cordinate array of first line
        :param y1: Y cordinate array of first line
        :param x2: X cordinate array of second line
        :param y2: Y cordinate array of second line
        '''
        # find the slope of line1
        line1_m1, i1, r1, p1, std1 = stats.linregress(x1, y1)

        # find the slope of line2
        line2_m2, i2, r2, p2, std2 = stats.linregress(x2, y2)

        # return angle
        return np.abs(np.rad2deg(np.arctan((line1_m1-line2_m2) /
                                           (1 + (line1_m1 * line2_m2)))))

    def airfoilplot(self, report_dir=None, show_plot=False):
        '''
        This method read the length of the output file object and set them to
        input plotting loop. It generates the following figures,
        figure 1: subplot of airfoil X,Y co-ordinates

        :param report_dir: base directory to which plot images to be written
        :param show_plot: a flag to display plot if needed.
        '''
        if report_dir is not None:
            plot_path = os.path.join(report_dir, "airfoil.jpg")

        # setting up canvas
        rcp.update({'font.size': 10})
        rcp['legend.loc'] = 'best'

        # set figure
        fig1 = plt.figure(figsize=(16.0, 9.0))
        fig1.canvas.set_window_title("Airfoil Profile Plot")
        # Airfoil X and Y
        ax1_fig1 = fig1.add_subplot(111, aspect='equal')
        ax1_fig1.grid(True)

        colour_mmap = ['#292421']

        # prep x and y
        data = np.vstack((self.airfoil_table, self.airfoil_table[0, :]))
                          ax1_fig1.plot(data[:, 0], data[:, 1], 
                          c=colour_mmap[0],
                          label=os.path.basename(self.airfoil_file),
                          marker='+')

        # ask matplotlib for the plotted objects and their labels
        lines, labels = ax1_fig1.get_legend_handles_labels()
        ax1_fig1.legend(lines, labels, loc=0)
        fig1.tight_layout()

        if report_dir is not None:
            try:
                fig1.savefig(plot_path)
            except Exception, err:
                self._myprint(str(err))
        if show_plot:
            try:
                fig1.show()
            except Exception, err:
                self._myprint(str(err))

    def catch_scons(self, target, source, env):
        '''
        This method catches the scons argument and executes the save plot
        options.

        :param target: A list of targets scons is trying to build
        :param source: A list of source files scons is using to build targets
        :param env: A list of all environmental varibles from scons.
        '''
        input_file = env['airfoil_profile_file']
        self.load_airfoilfile(input_file)
        report_dir = env['foil_report_dir']
        self.airfoilplot(report_dir)

    def _myprint(self, printstrg):
        '''
        This method is helper to log if a logger object is provide else prints
        to terminal.

        :param printstrg: string to be logged or print to terminal
        '''
        try:
            self.foil_io_logger.info(printstrg)
        except:
            print "AIRFOIL.IO : " + printstrg

    def get_foil_prop(self):
        '''
        This is a method to wrap all property method of the airfoils
        '''
        # calculate trailing Edge thick
        last_idx = len(self.airfoil_table)
        self.trail_thick = np.power((np.sum(np.power((self.airfoil_table
                                                      [last_idx - 1] -
                                                      self.airfoil_table[0]),
                                                      2))), 0.5)
        try:
            self.find_trail_angle()
        except:
            self._myprint("Error calculating trailing edge angle")
        try:
            self.thick_find()
        except:
            self._myprint("Error calculating thickness properties")
        try:
            self.find_pre_side_infl()
        except:
            self._myprint("Error calculating inflection point")
        try:
            self.find_LE_radius()
        except:
            self._myprint("Error calculating LE radius")
        try:
            self.find_airfoil_area()
        except:
            self._myprint("Error calculating airfoil area")
        try:
            self.find_camber()
        except:
            self._myprint("Error calculating camber")

if __name__ == "__main__":
    pass
