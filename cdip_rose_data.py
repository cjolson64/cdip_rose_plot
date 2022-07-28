import netCDF4
import numpy as np
import datetime


class cdip_rose_data:


   def __init__(self, station_id, start_date, delta_days, *kwargs):

      # Set up some defaults for the bin sizes
      self.radial_bin_count = 16
      self.height_bin_count = 10
      self.period_bin_count = 10

      self.height_max = 5 # meters
      self.period_max = 20 # seconds

      self.metric = True # change to False to convert to feet

      if (self.metric):
         self.length_convert_factor = 1.0
      else:
         self.length_convert_factor = 3.28084


      # Validate infput formats
      # start_date should be mm/dd/yyyy and delta_days should be int

      url = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/' + station_id + '/' + station_id + '_historic.nc'
      self.nc = netCDF4.Dataset(url)

      self.ncTime = self.nc.variables['waveTime'][:]
      self.Fq = self.nc.variables['waveFrequency'] # Assign variable name - Wave Frequency
      self.Tp = self.nc.variables['waveTp'] # Assign variable name - Peak Wave Period
      self.Dp = self.nc.variables['waveDp'] # Assign variable name - Peak Wave Direction
      self.Hs = self.nc.variables['waveHs'] # Assign variable name - Significant Wave Height

      start_unix = int(datetime.datetime.strptime(start_date, '%m/%d/%Y').timestamp())
      end_unix = start_unix + (delta_days * 86400)

      self.start_index = self.find_index_nearest_time(start_unix)
      self.end_index = self.find_index_nearest_time(end_unix)


      self._set_radial_bins()
      self._set_height_bins()
      self._set_period_bins()
      

   def _validate_input_date(start_date):

      return True


   def _set_radial_bins(self):

      self.radial_bin_width = 360/self.radial_bin_count
      self.radial_start = self.radial_bin_width/2
      self.radial_start_angles = np.arange(self.radial_start, 360, self.radial_bin_width)


   def _set_height_bins(self):

      self.height_bin_width = self.height_max/self.height_bin_count
      self.height_bin_starts = np.arange(0, self.height_max, self.height_bin_width)

   def _set_period_bins(self):

      self.period_bin_width = self.period_max/self.period_bin_count
      self.period_bin_starts = np.arange(0, self.period_max, self.period_bin_width)

   def find_index_nearest_time(self, unixtime):

      #nearest = find_nearest(self.ncTime, unixtime)  # Find the closest unix timestamp
      #near_index = np.where(self.ncTime==nearest)[0][0]  # Grab the index number of found date
      near_index = (np.abs(self.ncTime - unixtime)).argmin()

      return near_index
      

   def find_radial_bin_number(self, angle):

      return int((angle - self.radial_start) / self.radial_bin_width) - 1


   def find_height_bin_number(self, height):

      return int(height/self.height_bin_width) - 1


   def get_wave_height_rose_data(self):

      wave_directions = np.copy(self.Dp[self.start_index:self.end_index])
      wave_heights = np.copy(self.Hs[self.start_index:self.end_index])
      data = np.array([wave_directions, wave_heights]) 
      rose_counts = np.zeros([self.radial_bin_count, self.height_bin_count], 'int')

      for row in data.T:
         radial_bin = self.find_radial_bin_number(row[0])
         height_bin = self.find_height_bin_number(row[1])

         rose_counts[radial_bin, height_bin] += 1
          
      return rose_counts/np.sum(rose_counts)


   def get_wave_period_rose_data(self):

      wave_directions = np.copy(self.Dp[self.start_index:self.end_index])
      wave_heights = np.copy(self.Hs[self.start_index:self.end_index])
      data = np.array([wave_directions, wave_heights]) 
      rose_counts = np.zeros([self.radial_bin_count, self.height_bin_count], 'int')

      for row in data.T:
         radial_bin = self.find_radial_bin_number(row[0])
         height_bin = self.find_height_bin_number(row[1])

         rose_counts[radial_bin, height_bin] += 1
          
      return rose_counts/np.sum(rose_counts)




if __name__ == "__main__":

   station_id = '198p1'
   start_date = '11/01/2013'
   delta_days = 30

   rose_data = cdip_rose_data(station_id, start_date, delta_days)
   result = rose_data.get_wave_height_rose_data()

   print(result)


