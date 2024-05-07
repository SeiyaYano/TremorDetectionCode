import argparse
import pandas as pd
import os
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
import datetime
from obspy import *
import const
import modules
import numpy as np
from scipy import integrate, interpolate


def features_calculator(i,df):

    ymdhms = datetime.datetime(df["year"][i],df["month"][i],df["day"][i],df["hour"][i],df["min"][i],df["sec"][i],df["microsec"][i])
    ymd = str(ymdhms.year) + str(ymdhms.month).zfill(2) + str(ymdhms.day).zfill(2)
    year_month_day = datetime.datetime(ymdhms.year,ymdhms.month,ymdhms.day)
    p_utc = UTCDateTime(df["ymdhms"][i])

    hypo_pos = [df["lat"][i],df["lon"][i],df["dep"][i]]

    ch_list = eval(df["channel path"][i])
    ch_pos_list = eval(df["channel pos"][i])
    ch_list_in = []
    for _ in range(len(ch_pos_list)):
        ch_pos = ch_pos_list[_]
        res = modules.latlon2xy([hypo_pos[0],hypo_pos[1]],ch_pos)
        dis = res["distance"]
        if dis <= const.epi_dis:
            ch_list_in.append(ch_list[_])
    
    dt_list = []
    s_raw,s_band = Stream(),Stream()
    
    if ymdhms.hour == 23 and ymdhms.minute >= 60-const.timewindow_length*3/60:
        for j in range(len(ch_list_in)):
            if ymd[2:8] in str(ch_list_in[j]):
                ch_path_1 = ch_list_in[j]
                _ = year_month_day + datetime.timedelta(days=1)
                _d = str(_.year)[2:4] + str(_.month).zfill(2) + str(_.day).zfill(2)
                ch_path_2 = (ch_path_1.replace(str(ymdhms.year),str(_.year))).replace(ymd[2:8],_d)

            elif ymd[2:8] not in str(ch_list_in[j]):
                ch_path_1 = ch_list_in[j]
                _ = year_month_day + datetime.timedelta(days=1)
                _d = str(_.year)[2:4] + str(_.month).zfill(2) + str(_.day).zfill(2)
                ch_path_2 = (ch_path_1.replace(str(_.year),str(ymdhms.year))).replace(_d,ymd[2:8])

            s = None
            try:
                s = read(ch_path_1)
                s += read(ch_path_2)
                s.merge(method=1)
            except:
                f = open("errors2.txt","a")
                f.write('sac_read_error : %s - %s or %s\n' %(str(ymdhms),ch_path_1,ch_path_2))
                f.close()
                s = None
            
            if s is not None:
                st_la = s[0].stats.sac.stla
                st_lo = s[0].stats.sac.stlo
                ch_pos = [st_la,st_lo,0]
                dt = modules.japan_land_travel_time(ch_pos,hypo_pos,"S")
                dt_list.append(float(dt))
                t_s = p_utc + dt - df["dur pre"][i] - df["dur"][i]
                t_e = p_utc + dt + df["dur"][i]

                s_raw += s.slice(t_s,t_e)
                s.filter("bandpass",freqmin=const.low1,freqmax=const.low2,corners=const.corner,zerophase=True)
                s_band += s.slice(t_s,t_e)

    elif ymdhms.hour == 0 and ymdhms.minute <= const.timewindow_length*3/60:
        for j in range(len(ch_list_in)):
            ch_path_1 = ch_list_in[j]
            _ = year_month_day + datetime.timedelta(days=-1)
            _d = str(_.year)[2:4] + str(_.month).zfill(2) + str(_.day).zfill(2)
            ch_path_2 = (ch_path_1.replace(str(ymdhms.year),str(_.year))).replace(ymd[2:8],_d)
            
            s = None
            try:
                s = read(ch_path_1)
                s += read(ch_path_2)
                s.merge(method=1)
            except:
                f = open("errors2.txt","a")
                f.write('sac_read_error : %s - %s or %s\n' %(str(ymdhms),ch_path_1,ch_path_2))
                f.close()
                s = None

            if s is not None:
                st_la = s[0].stats.sac.stla
                st_lo = s[0].stats.sac.stlo
                ch_pos = [st_la,st_lo,0]
                dt = modules.japan_land_travel_time(ch_pos,hypo_pos,"S")
                dt_list.append(float(dt))
                t_s = p_utc + dt - df["dur pre"][i] - df["dur"][i]
                t_e = p_utc + dt + df["dur"][i]

                s_raw += s.slice(t_s,t_e)
                s.filter("bandpass",freqmin=const.low1,freqmax=const.low2,corners=const.corner,zerophase=True)
                s_band += s.slice(t_s,t_e)
    
    else:
        for j in range(len(ch_list_in)):
            s = None
            try:
                s = read(ch_list_in[j],starttime=p_utc-const.timewindow_length*3,endtime=p_utc+const.timewindow_length*3)
            except:
                f = open("errors2.txt","a")
                f.write('sac_read_error : %s - %s\n' %(str(ymdhms),ch_list[j]))
                f.close()
                s = None

            if s is not None:
                st_la = s[0].stats.sac.stla
                st_lo = s[0].stats.sac.stlo
                ch_pos = [st_la,st_lo,0]
                dt = modules.japan_land_travel_time(ch_pos,hypo_pos,"S")
                dt_list.append(float(dt))
                t_s = p_utc + dt - df["dur pre"][i] - df["dur"][i]
                t_e = p_utc + dt + df["dur"][i]

                s_raw += s.slice(t_s,t_e)
                s.filter("bandpass",freqmin=const.low1,freqmax=const.low2,corners=const.corner,zerophase=True)
                s_band += s.slice(t_s,t_e)

    N = df["dur"][i] * const.sampling_rate
    ch_pos_list_2 = []
    dur2s, cum_list = [], []
    Amp_e_list, Amp_n_list = [], []
    Amp_sn_list = []

    for j in range(len(s_raw)):
        s_utc = s_raw[j].times("utcdatetime")[0]
        sta_j,component_j = s_raw[j].stats.station,s_raw[j].stats.channel
        st_la_j = s_raw[j].stats.sac.stla
        st_lo_j = s_raw[j].stats.sac.stlo
        epi_dis_j = modules.latlon2xy([hypo_pos[0],hypo_pos[1]],[st_la_j,st_lo_j])["distance"]
        hypo_dis_j = (epi_dis_j**2 + df["dep"][i]**2)**0.5

        s_event = s_raw[j].slice(s_utc+df["dur"][i],s_utc+df["dur"][i]*2)
        s_event.detrend(type="linear")
        s_event.taper(0.1,type="cosine")
        y_fft_e = np.fft.fft(s_event.data)
        Amp_e = abs(y_fft_e/(N/2))[1:int(N/2)]

        freq_e = np.fft.fftfreq(N,d=1/const.sampling_rate)
        freq = freq_e[1:int(N/2)]

        s_noise = s_raw[j].slice(s_utc,s_utc+df["dur"][i])
        s_noise.detrend(type="linear")
        s_noise.taper(0.1,type="cosine")
        y_fft_n = np.fft.fft(s_noise.data)
        Amp_n = abs(y_fft_n/(N/2))[1:int(N/2)]

        Amp_sn_raw = np.array((Amp_e)**2 - (Amp_n)**2)
        Amp_sn_raw[Amp_sn_raw<0] = 0
        Amp_sn_raw = (Amp_sn_raw)**0.5
        SN = Amp_sn_raw/Amp_n 
        Amp_sn = Amp_sn_raw * hypo_dis_j

        Amp_e_list.append(Amp_e)
        Amp_n_list.append(Amp_n)
        Amp_sn[SN<const.sn_ratio] = None
        Amp_sn_list.append(Amp_sn)

        for k in range(len(freq)-1):
            if freq[k] <= const.low1 < freq[k+1]:
                if const.low1 in freq:
                    k_low1 = k
                else:
                    k_low1 = k+1
            if freq[k] <= const.low2 < freq[k+1]:
                k_low2 = k
            if freq[k] <= const.high1 < freq[k+1]:
                if const.high1 in freq:
                    k_high1 = k
                else:
                    k_high1 = k+1
            if freq[k] <= const.high2 < freq[k+1]:
                k_high2 = k

        d_low = pd.DataFrame(Amp_sn[k_low1:k_low2])
        d_high = pd.DataFrame(Amp_sn[k_high1:k_high2])

        if d_low.mean()[0]>0 or d_high.mean()[0]>0:
            ch_pos_list_2.append([st_la_j,st_lo_j])
            
            s_dur2 = s_band[j].slice(s_utc+df["dur pre"][i], s_utc+df["dur pre"][i]+df["dur"][i]*2)
            t = s_dur2.times()
            s_dur2 = s_dur2.detrend(type="linear")
            d = np.array(s_dur2)
            env = d**2
            int_cumulative_trapezoid = integrate.cumulative_trapezoid(env, t, initial=0)
            scaled_cum = int_cumulative_trapezoid/int_cumulative_trapezoid[-1]
#           cum_list.append(scaled_cum)
            dur2 = np.inf
            for x in range(0,501):
                x = x/1000
                y = x + 0.5
                cx = modules.inter(scaled_cum,t,"linear",x)
                cy = modules.inter(scaled_cum,t,"linear",y)
                if cy-cx < dur2:
                    dur2 = cy-cx
            dur2s.append(float(dur2)) 

    if len(dur2s) == 0:
        energy_dur, dur2s = -1,-1
        f_low_e, f_high_e, f_low_n, f_high_n, f_low_sn, f_high_sn = -1, -1, -1, -1, -1, -1
        f_ratio = -1
        ch_pos_list_2 = -1
    else:
        Amp_e_all, Amp_n_all, Amp_sn_all = pd.DataFrame(Amp_e_list), pd.DataFrame(Amp_n_list), pd.DataFrame(Amp_sn_list)
        Amp_e_median, Amp_n_median, Amp_sn_median = Amp_e_all.median(), Amp_n_all.median(), Amp_sn_all.median()
        f_low_e, f_high_e = float(Amp_e_median[k_low1:k_low2].mean()),float(Amp_e_median[k_high1:k_high2].mean())
        f_low_n, f_high_n = float(Amp_n_median[k_low1:k_low2].mean()),float(Amp_n_median[k_high1:k_high2].mean())
        f_low_sn, f_high_sn = float(Amp_sn_median[k_low1:k_low2].mean()),float(Amp_sn_median[k_high1:k_high2].mean())

        if f_high_sn > 0 and f_low_sn > 0:
            f_ratio = float(f_high_sn/f_low_sn)
            energy_dur = float(min(dur2s))
        else:
            energy_dur, dur2s = -1,-1
            f_low_e, f_high_e, f_low_n, f_high_n, f_low_sn, f_high_sn = -1, -1, -1, -1, -1, -1
            f_ratio = -1
            ch_pos_list_2 = -1

    a = [str(ymdhms),ymdhms.year,ymdhms.month,ymdhms.day,ymdhms.hour,ymdhms.minute,ymdhms.second,ymdhms.microsecond,df["doy"][i],
        df["lat"][i],df["lon"][i],df["mag"][i],df["Es"][i],df["dep"][i],df["dur"][i],df["dur pre"][i],df["dur pos"][i],
        df["channel path"][i],df["channel pos"][i],df["first az"][i],df["second az"][i],
        df["lat error"][i],df["lon error"][i],df["dep error"][i],df["cc"][i],
        energy_dur,dur2s,
        f_low_e,f_high_e,f_low_n,f_high_n,f_low_sn,f_high_sn,f_ratio,
        ch_pos_list_2]

    return a


def make_catalog(data_list):
    df = pd.DataFrame(data_list,columns=["ymdhms","year","month","day","hour","min","sec","microsec","doy",
                                        "lat","lon","mag","Es","dep","dur","dur pre","dur pos",
                                        "channel path","channel pos","first az","second az",
                                        "lat error","lon error","dep error","cc",
                                        "energy dur","energy dur list",
                                        "f low e","f high e","f low n","f high n","f low sn","f high sn","f ratio",
                                        "channel pos 2"])
    df = df[df["energy dur"]!=-1].copy()
    return df


def event_features(args):
    
    df = pd.read_csv(args.savepath + "/" + args.year + "/csv_files/events.csv")

    if os.path.isfile(args.savepath + "/" + args.year + "/csv_files/features.csv"):
        df2 = pd.read_csv(args.savepath + "/" + args.year + "/csv_files/features.csv")
        i_start = list(df["ymdhms"]).index(str(list(df2["ymdhms"])[-1])) + 1
    else:
        df2 = pd.DataFrame()
        i_start = 0

    print("event features calculation from ",i_start)
    df = df[i_start:].copy()
    df.reset_index(drop=True, inplace=True)
    with tqdm(total=len(df)) as progress:
        with ProcessPoolExecutor(max_workers=os.cpu_count()//2) as executor:
            futures = []
            for i in range(len(df)):
                future = executor.submit(features_calculator,i,df)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)
                
                if i != 0 and i % 1000 == 0:
                    data_list = [f.result() for f in futures]
                    df1 = make_catalog(data_list)
                    catalog = pd.concat([df1,df2])
                    catalog.sort_values(by="ymdhms", ascending=True, inplace=True)
                    catalog.reset_index(drop=True, inplace=True)
                    catalog.to_csv(args.savepath + "/" + args.year + "/csv_files/features.csv", index=False)
                
    data_list = [f.result() for f in futures]
    df1 = make_catalog(data_list)
    catalog = pd.concat([df1,df2])
    catalog.sort_values(by="ymdhms", ascending=True, inplace=True)
    catalog.reset_index(drop=True, inplace=True)
    catalog.to_csv(args.savepath + "/" + args.year + "/csv_files/features.csv", index=False)



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--year', help='year',required=True)
    parser.add_argument('--savepath', help='save path',required=True)
    args = parser.parse_args()

    event_features(args)