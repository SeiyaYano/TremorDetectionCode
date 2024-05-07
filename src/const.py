# Outlier control
c_lim = 0.6
c_lim_t = 0.4

# Error estimation
bootstrap_itr = 100
max_error = 2  # km

# Time window
timewindow_length = 300  # sec
timewindow_step = 150  # sec

# Source location
gridsearch_depth = 30  # km
min_source_depth = 0  # km
max_source_depth = 700  # km

# Search area
max_dist_station_source = 50  # km
max_source_depth_search = 500  # km
max_dist_pos_station = 200 #km
max_dist_from_grid = 200 #km

# Magnitude estimation
rho = 3000  # kg/m^3
beta = 2844  # m/s

# Detection criteria
min_pair_num = 16

# Station pair criteria
min_dist_station_pair = 0.3  # km
max_dist_station_pair = 100  # km

# Merge solutions
#merge_source_dist = 0.2  # degree

# Min duration
min_duration = 1  # s

#frequency band
freq_min = 2.0 #[Hz]
freq_max = 8.0 #[Hz]
low_pass = 0.2 #[Hz]
corner = 10.0 #[Hz]
low1 = 2.0 #[Hz]
low2 = 8.0 #[Hz]
high1 = 10.0 #[Hz]
high2 = 20.0 #[Hz]
sampling_rate = 100 #[sps]

#event quality control
azimuthal_gap = 180 #[km]
epi_dis = 50 #[km]
sn_ratio = 2.0 #[-]
probability = 0.99 #[-]

#figures
lat_range = [30.5,38]
lon_range = [129,138.5]

#clustering
min_s = 0 #[-]; specify positive integer or zero
fractal_dimenstion = 0
eta = 0
thre_T = 0
thre_R = 0