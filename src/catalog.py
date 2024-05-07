import modules


def make_catalog(time_window, num_station, time, cc, dur_pre, dur_pos, me, x, err, channel_path_list, channel_pos_list):

    station_pos_list = []
    for ch in channel_pos_list:
        ch_pos = [float(ch[0]),float(ch[1])]
        if ch_pos not in station_pos_list:
            station_pos_list.append(ch_pos)

    azi_list = []
    for st_pos in station_pos_list:
        res = modules.latlon2xy([float(x[0]),float(x[1])],st_pos)
        azi = res["azimuth"]
        azi_list.append(azi)
    azi_list = sorted(azi_list)
    
    gap = []
    for j in range(0,len(azi_list)):
        if j != len(azi_list) - 1:
            gap.append(azi_list[j+1]-azi_list[j])
        else:
            gap.append(azi_list[0]+360-azi_list[j])

    first_az = sorted(gap)[-1]
    second_az = sorted(gap)[-2]
    
    return {
        'time window': int(time_window),
        'number of pairs': int(num_station),
        'time': float(time),
        'cross correlation': float(cc),
        'duration pre': int(dur_pre),
        'duration pos': int(dur_pos),
        'duration': int(dur_pre+dur_pos),
        'me': float(me),
        'latitude': float(x[0]),
        'longitude': float(x[1]),
        'depth': float(x[2]),
        'latitude error': float(err[0]),
        'longitude error': float(err[1]),
        'depth error': float(err[2]),
        '1st az': float(first_az),
        '2nd az': float(second_az),
        'channel path': channel_path_list,
        'channel position': channel_pos_list
    }
