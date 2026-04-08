"""
Calculate stresses over a rectangular geographic region on a regular
lat-lon-time grid.

L{gridcalc} allows you to specify a grid of times and locations (using a
L{Grid} object) at which to perform stress calculations using the L{satstress}
module.  A L{GridCalc} object is created from a L{StressCalc} object and a
L{Grid} object, and can be saved to disk for visualization and sharing as a
U{Unidata NetCDF data cube <http://www.unidata.ucar.edu/software/netcdf>} (.nc)
file.

"""

#Imports and Constants

import re
import time
import numpy
import scipy
import satstress as ss
import physcon as pc
import netCDF4 as netCDF3  # modern replacement for netCDF3
from optparse import OptionParser
from satstress.satstress import nvf2dict

__NETCDF_CONVENTIONS__ = "None"

#Main Execution
def main():
    usage = "usage: %prog [options] satfile gridfile outfile"
    description = __doc__
    op = OptionParser(usage)
    (options, args) = op.parse_args()

    #do a little error checking
    if len(args) != 3:
        op.error("incorrect number of arguments")

    with open(args[0], 'r') as satfile:
        the_sat = ss.Satellite(satfile)

    with open(args[1], 'r') as gridfile:
        the_grid = Grid(gridfile, satellite=the_sat)

    # Set up the stress model and calculate grid-based stress tensors
    the_stresscalc = ss.StressCalc([ss.NSR(the_sat), ss.Diurnal(the_sat)])
    the_gridcalc = GridCalc(the_grid, the_stresscalc)
    the_gridcalc.write_netcdf(args[2])

#Grid Definition

class Grid(object): # {{{
    """
    A container class defining the temporal and geographic range and resolution
    of the gridded stress calculation.

    The parameters defining the calculation grid are read in from a name value
    file, parsed into a Python dictionary using L{nvf2dict}, and used
    to set the data attributes of the L{Grid} object.

    The geographic extent and resolution of the calculation is specified by
    minimum and maximum values for latitude and longitude, as well as the
    number of regularly spaced latitude and longitude values to calcualte in
    total (including the minimum and maximum values).  Similarly, time or
    orbital position are specified by an initial value, a final value, and the
    total number of values the time dimension should take on.  In the input
    name-value file, the names used to specify these grid parameters are:

    LAT_MIN, LAT_MAX, LAT_NUM
    LON_MIN, LON_MAX, LON_NUM

    and only one of the following two temporal parameters:

    TIME_MIN,  TIME_MAX,  TIME_NUM
    ORBIT_MIN, ORBIT_MAX, ORBIT_NUM

    with latitude, longitude, and orbital position being specified in degrees,
    and time being specified in seconds.  Both orbital position and time are
    assumed to be zero at periapse.

    The variability of the NSR stresses are similarly explored by specifying 
    the range and resolution of the NSR_PERIOD using:

    NSR_PERIOD_MIN, NSR_PERIOD_MAX, NSR_PERIOD_NUM

    @ivar grid_id: A string identifying the grid
    @type grid_id: str

    @ivar lat_min: Southern bound, degrees (north positive).
    @type lat_min: float
    @ivar lat_max: Northern bound, degrees (north positive).
    @type lat_max: float
    @ivar lat_num: Number of latitude values in the grid.
    @type lat_num: float

    @ivar lon_min: Western bound, degrees (east positive).
    @type lon_min: float
    @ivar lon_max: Eastern bound, degrees (east positive).
    @type lon_max: float
    @ivar lon_num: Number of longitude values in the grid.
    @type lon_num: float

    @ivar time_min: Initial time at which calculation begins (0 = periapse).
    @type time_min: float
    @ivar time_max: Final time at which calculation ends.
    @type time_max: float
    @ivar time_num: Number of timesteps to calculate.
    @type time_num float

    @ivar orbit_min: Initial orbital position in degrees (0 = periapse)
    @type orbit_min: float
    @ivar orbit_max: Final orbital position in degrees (0 = periapse)
    @type orbit_max: float
    @ivar orbit_num: Number of orbital timesteps to calculate.
    @type orbit_num: float

    @ivar satellite: the satellite whose orbital parameters should be used in
    converting between orbital position and time (if necessary)
    @type satellite: L{satstress.Satellite}

    """

    def __init__(self, gridfile=None, satellite=None):
        """Initialize the Grid object from a gridfile.

        @param gridfile: a name value file specifying a calculation grid
        @type gridfile: file
        @param satellite: the satellite whose orbital parameters should be used
        in converting between orbital position and time (if necessary)
        @type satellite: L{satstress.Satellite}

        @raise MissingDimensionError: if the input gridfile does not specify
        the range and resolution of latitude, longitude, time/orbital position
        and NSR period values to do calculations for.

        """

        print(">>> Grid.__init__ called")
        assert gridfile is not None and satellite is not None
        gridParams = nvf2dict(gridfile, comment='#')

        self.satellite = satellite
        self.grid_id = gridParams['GRID_ID']

        try:
            self.lat_min = float(gridParams['LAT_MIN'])
            self.lat_max = float(gridParams['LAT_MAX'])
            self.lat_num = float(gridParams['LAT_NUM'])
        except KeyError:
            raise MissingDimensionError(gridfile, 'latitude')

        try:
            self.lon_min = float(gridParams['LON_MIN'])
            self.lon_max = float(gridParams['LON_MAX'])
            self.lon_num = float(gridParams['LON_NUM'])
        except KeyError:
            raise MissingDimensionError(gridfile, 'longitude')

        try:
            self.time_min = float(gridParams['TIME_MIN'])
            self.time_max = float(gridParams['TIME_MAX'])
            self.time_num = float(gridParams['TIME_NUM'])
            self.orbit_min = None
            self.orbit_max = None
            self.orbit_num = None
        except KeyError:
            try:
                self.orbit_min = float(gridParams['ORBIT_MIN'])
                self.orbit_max = float(gridParams['ORBIT_MAX'])
                self.orbit_num = float(gridParams['ORBIT_NUM'])
                self.time_min = satellite.orbit_period() * (self.orbit_min / 360.0)
                self.time_max = satellite.orbit_period() * (self.orbit_max / 360.0)
                self.time_num = satellite.orbit_period() * (self.orbit_num / 360.0)
            except KeyError:
                raise MissingDimensionError(gridfile, 'time/orbital position')

        try:
            self.nsr_period_min = float(gridParams['NSR_PERIOD_MIN'])
            self.nsr_period_max = float(gridParams['NSR_PERIOD_MAX'])
            self.nsr_period_num = float(gridParams['NSR_PERIOD_NUM'])
        except KeyError:
            raise MissingDimensionError(gridfile, 'NSR period')

    def __str__(self):
        """Output a grid definition file."""
        return f"<Grid {self.grid_id}: Lat {self.lat_min}-{self.lat_max}, Lon {self.lon_min}-{self.lon_max}>"

# GridCalc Engine

class GridCalc:
    """
    GridCalc coordinates the application of a stress model over a defined
    lat-lon-time grid, storing results and managing output.

        An object that performs a L{StressCalc} on a regularly spaced L{Grid}.

    A C{GridCalc} object takes a particular L{StressCalc} object and
    instantiates the calculation it embodies at each point in the regularly
    spaced grid specified by the associated L{Grid} object.

    """

    def __init__(self, grid=None, stresscalc=None):
        if grid is None or stresscalc is None:
            raise GridCalcInitError()
        self.grid = grid
        self.stresscalc = stresscalc

    def __str__(self):
        return str(self.grid.satellite) + str(self.grid)

    def write_netcdf(self, outfile):
        """
        Output a netCDF file containing the results of the calculation
        specified by the GridCalc object.

        Each stress field encapsulated in the GridCalc object will be output
        within the netCDF file as three data fields, one for each of the stress
        tensor components L{Ttt}_NAME, L{Tpt}_NAME, L{Tpp}_NAME, where NAME is
        the name of the L{StressDef} object (e.g. L{Diurnal} or L{NSR}).

        Writing out the calculation results causes the calculation to take
        place.  No mechanism for performing the calculation and retaining it
        in memory for manipulation is currently provided.

        """

        # Create a netCDF file object to stick the calculation results in:
        nc_out = netCDF3.Dataset(outfile, 'w')

        # Set metadata fields of nc_out appropriate to the calculation at hand.

        nc_out.description = "satstress calculation on a regular grid.  All parameter units are SI (meters-kilograms-seconds)"
        nc_out.history     = """Created: %s using the satstress python package: http://code.google.com/p/satstress""" % ( time.ctime(time.time()))
        nc_out.Conventions = __NETCDF_CONVENTIONS__

                # Independent parameters
        nc_out.grid_id = self.grid.grid_id
        nc_out.system_id = self.grid.satellite.system_id
        nc_out.planet_mass = self.grid.satellite.planet_mass
        nc_out.orbit_eccentricity = self.grid.satellite.orbit_eccentricity
        nc_out.orbit_semimajor_axis = self.grid.satellite.orbit_semimajor_axis

        for i in range(4):
            layer = self.grid.satellite.layers[i]
            setattr(nc_out, f"layer_id_{i}", layer.layer_id)
            setattr(nc_out, f"density_{i}", layer.density)
            setattr(nc_out, f"lame_mu_{i}", layer.lame_mu)
            setattr(nc_out, f"lame_lambda_{i}", layer.lame_lambda)
            setattr(nc_out, f"thickness_{i}", layer.thickness)
            setattr(nc_out, f"viscosity_{i}", layer.viscosity)

        # Derived properties
        nc_out.satellite_radius = self.grid.satellite.radius()
        nc_out.satellite_mass = self.grid.satellite.mass()
        nc_out.satellite_density = self.grid.satellite.density()
        nc_out.satellite_surface_gravity = self.grid.satellite.surface_gravity()
        nc_out.satellite_orbit_period = self.grid.satellite.orbit_period()

        # Placeholder metadata for SSWeb interface
        nc_out.ssweb_run_id = ""
        nc_out.ssweb_username = ""
        nc_out.ssweb_ip_address = ""

        # Define coordinate axes
        nc_out.createDimension('latitude', int(self.grid.lat_num))
        lats = nc_out.createVariable('latitude', 'f4', ('latitude',))
        lats.units = "degrees_north"
        lats.long_name = "latitude"
        lats[:] = numpy.linspace(self.grid.lat_min, self.grid.lat_max, int(self.grid.lat_num))

        nc_out.createDimension('longitude', int(self.grid.lon_num))
        lons = nc_out.createVariable('longitude', 'f4', ('longitude',))
        lons.units = "degrees_east"
        lons.long_name = "longitude"
        lons[:] = numpy.linspace(self.grid.lon_min, self.grid.lon_max, int(self.grid.lon_num))

        nc_out.createDimension('nsr_period', int(self.grid.nsr_period_num))
        nsr_periods = nc_out.createVariable('nsr_period', 'f4', ('nsr_period',))
        nsr_periods.units = "seconds"
        nsr_periods.long_name = "NSR period"
        nsr_periods[:] = numpy.logspace(numpy.log10(self.grid.nsr_period_min), numpy.log10(self.grid.nsr_period_max), int(self.grid.nsr_period_num))

        if self.grid.orbit_min is None:
            nc_out.createDimension('time', int(self.grid.time_num))
            times = nc_out.createVariable('time', 'f4', ('time',))
            times.units = "seconds"
            times.long_name = "time after periapse"
            times[:] = numpy.linspace(self.grid.time_min, self.grid.time_max, int(self.grid.time_num))
        else:
            nc_out.createDimension('time', int(self.grid.orbit_num))
            times = nc_out.createVariable('time', 'f4', ('time',))
            times.units = "degrees"
            times.long_name = "degrees after periapse"
            times[:] = numpy.linspace(self.grid.orbit_min, self.grid.orbit_max, int(self.grid.orbit_num))
        # Independent parameters
        nc_out.grid_id = self.grid.grid_id
        nc_out.system_id = self.grid.satellite.system_id
        nc_out.planet_mass = self.grid.satellite.planet_mass
        nc_out.orbit_eccentricity = self.grid.satellite.orbit_eccentricity
        nc_out.orbit_semimajor_axis = self.grid.satellite.orbit_semimajor_axis

        for i in range(4):
            layer = self.grid.satellite.layers[i]
            setattr(nc_out, f"layer_id_{i}", layer.layer_id)
            setattr(nc_out, f"density_{i}", layer.density)
            setattr(nc_out, f"lame_mu_{i}", layer.lame_mu)
            setattr(nc_out, f"lame_lambda_{i}", layer.lame_lambda)
            setattr(nc_out, f"thickness_{i}", layer.thickness)
            setattr(nc_out, f"viscosity_{i}", layer.viscosity)

        # Derived properties
        nc_out.satellite_radius = self.grid.satellite.radius()
        nc_out.satellite_mass = self.grid.satellite.mass()
        nc_out.satellite_density = self.grid.satellite.density()
        nc_out.satellite_surface_gravity = self.grid.satellite.surface_gravity()
        nc_out.satellite_orbit_period = self.grid.satellite.orbit_period()

        # Placeholder metadata for SSWeb interface
        nc_out.ssweb_run_id = ""
        nc_out.ssweb_username = ""
        nc_out.ssweb_ip_address = ""

        # Define coordinate axes
        nc_out.createDimension('latitude', int(self.grid.lat_num))
        lats = nc_out.createVariable('latitude', 'f4', ('latitude',))
        lats.units = "degrees_north"
        lats.long_name = "latitude"
        lats[:] = numpy.linspace(self.grid.lat_min, self.grid.lat_max, int(self.grid.lat_num))

        nc_out.createDimension('longitude', int(self.grid.lon_num))
        lons = nc_out.createVariable('longitude', 'f4', ('longitude',))
        lons.units = "degrees_east"
        lons.long_name = "longitude"
        lons[:] = numpy.linspace(self.grid.lon_min, self.grid.lon_max, int(self.grid.lon_num))

        nc_out.createDimension('nsr_period', int(self.grid.nsr_period_num))
        nsr_periods = nc_out.createVariable('nsr_period', 'f4', ('nsr_period',))
        nsr_periods.units = "seconds"
        nsr_periods.long_name = "NSR period"
        nsr_periods[:] = numpy.logspace(numpy.log10(self.grid.nsr_period_min), numpy.log10(self.grid.nsr_period_max), int(self.grid.nsr_period_num))

        if self.grid.orbit_min is None:
            nc_out.createDimension('time', int(self.grid.time_num))
            times = nc_out.createVariable('time', 'f4', ('time',))
            times.units = "seconds"
            times.long_name = "time after periapse"
            times[:] = numpy.linspace(self.grid.time_min, self.grid.time_max, int(self.grid.time_num))
        else:
            nc_out.createDimension('time', int(self.grid.orbit_num))
            times = nc_out.createVariable('time', 'f4', ('time',))
            times.units = "degrees"
            times.long_name = "degrees after periapse"
            times[:] = numpy.linspace(self.grid.orbit_min, self.grid.orbit_max, int(self.grid.orbit_num))
            
        # At this point, we should have all the netCDF dimensions and their
        # corresponding coordinate variables created (latitutde, longitude,
        # time/orbit, nsr_period), but we still haven't created the data
        # variables, which will ultimately hold the results of our stress
        # calculation, and which depend on the aforedefined dimensions


       # Create Diurnal stress tensor variables
        Ttt_Diurnal = nc_out.createVariable('Ttt_Diurnal', 'f4', ('time', 'latitude', 'longitude',))
        Ttt_Diurnal.units = "Pa"
        Ttt_Diurnal.long_name = "north-south component of Diurnal eccentricity stresses"

        Tpt_Diurnal = nc_out.createVariable('Tpt_Diurnal', 'f4', ('time', 'latitude', 'longitude',))
        Tpt_Diurnal.units = "Pa"
        Tpt_Diurnal.long_name = "shear component of Diurnal eccentricity stresses"

        Tpp_Diurnal = nc_out.createVariable('Tpp_Diurnal', 'f4', ('time', 'latitude', 'longitude',))
        Tpp_Diurnal.units = "Pa"
        Tpp_Diurnal.long_name = "east-west component of Diurnal eccentricity stresses"

        # Create NSR stress tensor variables
        Ttt_NSR = nc_out.createVariable('Ttt_NSR', 'f4', ('nsr_period', 'latitude', 'longitude',))
        Ttt_NSR.units = "Pa"
        Ttt_NSR.long_name = "north-south component of NSR stresses"

        Tpt_NSR = nc_out.createVariable('Tpt_NSR', 'f4', ('nsr_period', 'latitude', 'longitude',))
        Tpt_NSR.units = "Pa"
        Tpt_NSR.long_name = "shear component of NSR stresses"

        Tpp_NSR = nc_out.createVariable('Tpp_NSR', 'f4', ('nsr_period', 'latitude', 'longitude',))
        Tpp_NSR.units = "Pa"
        Tpp_NSR.long_name = "east-west component of NSR stresses"

        # Separate out Diurnal and NSR objects
        for stress in self.stresscalc.stresses:
            if stress.__name__ == 'Diurnal':
                diurnal_stress = ss.StressCalc([stress])
            if stress.__name__ == 'NSR':
                nsr_stress = ss.StressCalc([stress])

        # Loop through times and fill in Diurnal stress tensors
        for t in range(len(times[:])):
            if self.grid.orbit_min is None:
                time_sec = times[t]
            else:
                time_sec = diurnal_stress.stresses[0].satellite.orbit_period() * (times[t] / 360.0)

            print(f"Calculating Diurnal stresses at {times[t]} {times.long_name}")
            for lon in range(len(lons[:])):
                for lat in range(len(lats[:])):
                    Tau_D = diurnal_stress.tensor(theta=scipy.radians(90.0 - lats[lat]),
                                                  phi=scipy.radians(lons[lon]),
                                                  t=time_sec)
                    Ttt_Diurnal[t, lat, lon] = Tau_D[0, 0]
                    Tpt_Diurnal[t, lat, lon] = Tau_D[1, 0]
                    Tpp_Diurnal[t, lat, lon] = Tau_D[1, 1]

        nc_out.sync()

        # Zero eccentricity for NSR
        nsr_stress.stresses[0].satellite.orbit_eccentricity = 0.0

        for p_nsr in range(len(nsr_periods[:])):
            new_sat = nsr_stress.stresses[0].satellite
            new_sat.nsr_period = nsr_periods[p_nsr]
            nsr_stress = ss.StressCalc([ss.NSR(new_sat)])

            print(f"Calculating NSR stresses for Pnsr = {nsr_periods[p_nsr]} {nsr_periods.units}")
            for lon in range(len(lons[:])):
                for lat in range(len(lats[:])):
                    Tau_N = nsr_stress.tensor(theta=scipy.radians(90 - lats[lat]),
                                              phi=scipy.radians(lons[lon]),
                                              t=0)
                    Ttt_NSR[p_nsr, lat, lon] = Tau_N[0, 0]
                    Tpt_NSR[p_nsr, lat, lon] = Tau_N[1, 0]
                    Tpp_NSR[p_nsr, lat, lon] = Tau_N[1, 1]

        nc_out.sync()


        # Zero eccentricity for NSR
        nsr_stress.stresses[0].satellite.orbit_eccentricity = 0.0

        for p_nsr in range(len(nsr_periods[:])):
            new_sat = nsr_stress.stresses[0].satellite
            new_sat.nsr_period = nsr_periods[p_nsr]
            nsr_stress = ss.StressCalc([ss.NSR(new_sat)])

            print(f"Calculating NSR stresses for Pnsr = {nsr_periods[p_nsr]} {nsr_periods.units}")
            for lon in range(len(lons[:])):
                for lat in range(len(lats[:])):
                    Tau_N = nsr_stress.tensor(theta=scipy.radians(90 - lats[lat]),
                                              phi=scipy.radians(lons[lon]),
                                              t=0)
                    Ttt_NSR[p_nsr, lat, lon] = Tau_N[0, 0]
                    Tpt_NSR[p_nsr, lat, lon] = Tau_N[1, 0]
                    Tpp_NSR[p_nsr, lat, lon] = Tau_N[1, 1]

        nc_out.sync()

#Error Types

class Error(Exception):
    """Base class for gridcalc-specific exceptions."""
    pass

class GridParamError(Error):
    """Raised when there is an issue with grid parameter values."""
    pass

class GridCalcInitError(Error):
    """Raised when GridCalc cannot initialize properly."""
    def __str__(self):
        return "GridCalc could not be initialized due to missing grid or stresscalc."

class MissingDimensionError(GridParamError):
    """
    Raised when a required grid dimension (latitude, longitude, time, or NSR period)
    is missing from the grid input file.
    """
    def __init__(self, gridfile, missing_dim):
        self.gridfile = gridfile
        self.missing_dim = missing_dim

    def __str__(self):
        return (f"\nNo range of {self.missing_dim} values was specified in the calculation grid definition\n"
                f"from the file:\n\n{self.gridfile.name}\n\nEvery Grid must contain at least a single time value or orbital position.\n")

#Run Script Entry Point

if __name__ == "__main__":
    main()

# What do I actually want the plotting tool to do?
#
# Given a list of netCDF files containing variables Ttt_*, Tpt_*, Tpp_*
# representing the tensor components of the surface membrane stresses:

#  - allow the display of either the magnitude of the tensor component,
#  - or the magnitude and direction of the principal components
#  - of one, or the sum of several, of the stress tensors in the files.

# Assuming that all of the netCDF files given to the tool contain compatible
# grids (the same set of lat/lon points) the user should be able to choose
# within each file amongst the values available for that file's "record"
# dimension (e.g. the time/orbital position dimension for Diurnal)
#
# Ideally, the user would be able to use even grids that don't have the same
# geographic points, so long as they are close enough in density for an
# interpolated grid to be acceptably accurate, but that's gravy, and something
# that would be part of the plotting program, not gridcalc.
