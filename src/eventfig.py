import argparse
import pandas as pd
from obspy import *
import datetime
import modules
import const
import numpy as np
from scipy import integrate, interpolate
import subprocess
import matplotlib.pyplot as plt


def fig_save(f,args,name):
    if str(type(f)) == "<class 'pygmt.figure.Figure'>":
        if args.outfile == "png":
            f.savefig(args.savepath + "/eventfig/" + name + ".png")
        elif args.outfile == "eps":
            f.savefig(args.savepath + "/eventfig/eps/" + name + ".eps")
    elif str(type(f)) == "<class 'obspy.core.stream.Stream'>":
        if args.outfile == "png":
            f.plot(size=(800,600), outfile = args.savepath + "/eventfig/" + name + ".png")
        elif args.outfile == "eps":
            f.plot(size=(800,600), outfile = args.savepath + "/eventfig/eps/" + name + ".eps")
    elif str(type(f)) == "<class 'module'>":
        if args.outfile == "png":
            f.savefig(args.savepath + "/eventfig/" + name + ".png",bbox_inches="tight",pad_inches=0.05)
        elif args.outfile == "eps":
            f.savefig(args.savepath + "/eventfig/eps/" + name + ".eps",bbox_inches="tight",pad_inches=0.05)
    elif str(type(f)) == "<class 'PIL.Image.Image'>":
        if args.outfile == "png":
            f.save(args.savepath + "/eventfig/" + name + ".png")
        elif args.outfile == "eps":
            f.save(args.savepath + "/eventfig/eps/" + name + ".eps")


def event_fig(args,df):
    
    i = list(df["ymdhms"]).index(str(args.ymdhms))
    
    year,month,day = int(df["year"][i]),int(df["month"][i]),int(df["day"][i])
    hour,minu,sec,microsec = int(df["hour"][i]),int(df["min"][i]),int(df["sec"][i]),int(df["microsec"][i])
    p_utc = UTCDateTime("%s-%s-%sT%s:%s:%s.%sZ"
                        %(year,str(month).zfill(2),str(day).zfill(2),
                          str(hour).zfill(2),str(minu).zfill(2),str(sec).zfill(2),str(microsec).zfill(2)))

    ymdhms = datetime.datetime(year,month,day,hour,minu,sec,microsec)
    date = datetime.datetime(year=year,month=month,day=day)
    ymd = str(year) + str(month).zfill(2) + str(day).zfill(2)
    date_str = ymd[2:8]
    
    lat,lon,dep = float(df["lat"][i]),float(df["lon"][i]),float(df["dep"][i])
    dur,dur_pre,dur_pos = int(df["dur"][i]),int(df["dur pre"][i]),int(df["dur pos"][i])

    ch_list = eval(df["channel path"][i])
    ch_pos_list = eval(df["channel pos"][i])

    ch_list_in = []
    for _ in range(len(ch_pos_list)):
        ch_pos = ch_pos_list[_]
        res = modules.latlon2xy([lat,lon],ch_pos)
        dis = res["distance"]
        if dis <= const.epi_dis:
            ch_list_in.append(ch_list[_])
    
    map_la = [lat-2,lat+2]
    map_lo = [lon-3.5,lon+3.5]
    
    ##########
    dt_list = []
    ch_pos_list_2 = []
    
    s_raw,s_band = Stream(),Stream()

    if ymdhms.hour == 23 and ymdhms.minute >= 60-const.timewindow_length*3/60:
        for j in range(len(ch_list_in)):
            if ymd[2:8] in str(ch_list_in[j]):
                ch_path_1 = ch_list_in[j]
                _ = date + datetime.timedelta(days=1)
                _d = str(_.year)[2:4] + str(_.month).zfill(2) + str(_.day).zfill(2)
                ch_path_2 = (ch_path_1.replace(str(ymdhms.year),str(_.year))).replace(ymd[2:8],_d)

            elif ymd[2:8] not in str(ch_list_in[j]):
                ch_path_1 = ch_list_in[j]
                _ = date + datetime.timedelta(days=1)
                _d = str(_.year)[2:4] + str(_.month).zfill(2) + str(_.day).zfill(2)
                ch_path_2 = (ch_path_1.replace(str(_.year),str(ymdhms.year))).replace(_d,ymd[2:8])

            s = None
            try:
                s = read(ch_path_1)
                s += read(ch_path_2)
                s.merge(method=1)
            except:
                s = None
            
            if s is not None:
                st_la = s[0].stats.sac.stla
                st_lo = s[0].stats.sac.stlo
                ch_pos = [st_la,st_lo,0]
                dt = modules.japan_land_travel_time(ch_pos,[lat,lon,dep],"S")
                dt_list.append(float(dt))
                t_s = p_utc + dt - df["dur pre"][i] - df["dur"][i]
                t_e = p_utc + dt + df["dur"][i]

                s_raw += s.slice(t_s,t_e)
                s.filter("bandpass",freqmin=const.low1,freqmax=const.low2,corners=const.corner,zerophase=True)
                s_band += s.slice(t_s,t_e)

    elif ymdhms.hour == 0 and ymdhms.minute <= const.timewindow_length*3/60:
        for j in range(len(ch_list_in)):
            ch_path_1 = ch_list_in[j]
            _ = date + datetime.timedelta(days=-1)
            _d = str(_.year)[2:4] + str(_.month).zfill(2) + str(_.day).zfill(2)
            ch_path_2 = (ch_path_1.replace(str(ymdhms.year),str(_.year))).replace(ymd[2:8],_d)
            
            s = None
            try:
                s = read(ch_path_1)
                s += read(ch_path_2)
                s.merge(method=1)
            except:
                s = None

            if s is not None:
                st_la = s[0].stats.sac.stla
                st_lo = s[0].stats.sac.stlo
                ch_pos = [st_la,st_lo,0]
                dt = modules.japan_land_travel_time(ch_pos,[lat,lon,dep],"S")
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
                s = None

            if s is not None:
                st_la = s[0].stats.sac.stla
                st_lo = s[0].stats.sac.stlo
                ch_pos = [st_la,st_lo,0]
                dt = modules.japan_land_travel_time(ch_pos,[lat,lon,dep],"S")
                dt_list.append(float(dt))
                t_s = p_utc + dt - df["dur pre"][i] - df["dur"][i]
                t_e = p_utc + dt + df["dur"][i]

                s_raw += s.slice(t_s,t_e)
                s.filter("bandpass",freqmin=const.low1,freqmax=const.low2,corners=const.corner,zerophase=True)
                s_band += s.slice(t_s,t_e)

    fig_save(s_band,args,ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"waveform")

    N = dur * const.sampling_rate
    sta_chs = []
    dur2s, cum_list = [], []
    Amp_e_list, Amp_n_list = [], []
    Amp_sn_list = []
    ts = []
    
    for j in range(len(s_raw)):
        s_utc = s_raw[j].times("utcdatetime")[0]
        sta_j,component_j = s_raw[j].stats.station,s_raw[j].stats.channel
        st_la_j = s_raw[j].stats.sac.stla
        st_lo_j = s_raw[j].stats.sac.stlo
        epi_dis_j = modules.latlon2xy([lat,lon],[st_la_j,st_lo_j])["distance"]
        hypo_dis_j = (epi_dis_j**2 + dep**2)**0.5
        
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
            sta_chs.append("%s.%s" %(sta_j,component_j))
            s_dur2 = s_band[j].slice(s_utc+df["dur pre"][i], s_utc+df["dur pre"][i]+df["dur"][i]*2)
            t = s_dur2.times()
            ts.append(len(t))
            s_dur2 = s_dur2.detrend(type="linear")
            d = np.array(s_dur2)
            env = d**2
            int_cumulative_trapezoid = integrate.cumulative_trapezoid(env, t, initial=0)
            scaled_cum = int_cumulative_trapezoid/int_cumulative_trapezoid[-1]
            cum_list.append(scaled_cum)
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
        ch_pos_l = np.asarray(ch_pos_list)
        fig = modules.hypo_mono(ch_pos_l[:,0],ch_pos_l[:,1],c="darkgray",style="i0.15c",pen="black",lat_range=map_la,lon_range=map_lo)
        fig.plot(x=lon,y=lat,fill="red",style="a0.2c",pen="black")
        if args.text == "Y":
            text_list = [[map_la[1]-0.1,map_lo[0]+0.2,"%s (s) - %s - %s (s)" %(dur_pre,ymdhms,dur_pos)],
                         [map_la[0]+1.28,map_lo[1]-3,"depth %s (km)" %dep],
                         [map_la[0]+1.08,map_lo[1]-3,"envelope dur %s (s)" %dur],
                         [map_la[0]+0.88,map_lo[1]-3,"energy dur - (s)"],
                         [map_la[0]+0.68,map_lo[1]-3,"magnitude %s" %df["mag"][i]],
                         [map_la[0]+0.48,map_lo[1]-3,"1st azimuthal gap %s" %df["first az"][i]],
                         [map_la[0]+0.28,map_lo[1]-3,"2nd azimuthal gap %s" %df["second az"][i]],
                         [map_la[0]+0.08,map_lo[1]-3,"ratio -"]]
            fig = modules.map_text(text_list,font="8p",angle=0,justify="ML",fig=fig)
            fig_save(fig,args,ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"map")
            if args.outfile == "png":
                im1_path = args.savepath+"/eventfig/"+ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"map"+".png"
                im2_path = args.savepath+"/eventfig/"+ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"waveform"+".png"
                pic = modules.get_concat_h_resize(im1_path, im2_path)
                pic.save(args.savepath+"/eventfig/"+ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+".png")
    else:
        Amp_e_all, Amp_n_all, Amp_sn_all = pd.DataFrame(Amp_e_list), pd.DataFrame(Amp_n_list), pd.DataFrame(Amp_sn_list)
        Amp_e_median, Amp_n_median, Amp_sn_median = Amp_e_all.median(), Amp_n_all.median(), Amp_sn_all.median()
        f_low_e, f_high_e = float(Amp_e_median[k_low1:k_low2].mean()),float(Amp_e_median[k_high1:k_high2].mean())
        f_low_n, f_high_n = float(Amp_n_median[k_low1:k_low2].mean()),float(Amp_n_median[k_high1:k_high2].mean())
        f_low_sn, f_high_sn = float(Amp_sn_median[k_low1:k_low2].mean()),float(Amp_sn_median[k_high1:k_high2].mean())

        if f_high_sn > 0 and f_low_sn > 0:
            f_ratio = float(f_high_sn/f_low_sn)
            energy_dur = float(min(dur2s))
        
        #dur2
        t_len = min(ts)        
        s_c_n = dur2s.index(energy_dur)
        pfa = [None]*3
        for cum in cum_list:
            pfa = modules.simple_graph(t[0:t_len],cum[0:t_len],c="lightgray",pfa=pfa)
        pfa = modules.simple_graph(t[0:t_len],cum_list[s_c_n][0:t_len],xlabel="Time (s)",ylabel="Scaled cumulative energy",c="orangered",pfa=pfa)
        pfa[2].text(0,0.95,sta_chs[s_c_n])
        fig_save(pfa[0],args,ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"cum")

        #spectrum
        pfa = [None]*3
        for s in Amp_n_list:
            pfa = modules.simple_graph(freq,s,c="lightgray",lw=0.25,pfa=pfa)
        for s in Amp_e_list:
            pfa = modules.simple_graph(freq,s,c="skyblue",lw=0.25,pfa=pfa)
        pfa = modules.simple_graph(freq,Amp_n_median,c="black",pfa=pfa)
        pfa = modules.simple_graph(freq,Amp_e_median,c="royalblue",xlabel="Frequency (Hz)",ylabel="Amplitude (m/s)",xlog=True,ylog=True,pfa=pfa)
        fig_save(pfa[0],args,ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"spe")

        #spectrum (event - BG)
        pfa = [None]*3
        for s in Amp_sn_list:
            pfa = modules.simple_graph(freq,s,c="lightgray",lw=0.5,pfa=pfa)
        pfa = modules.simple_graph(freq,Amp_sn_median,c="black",xlabel="Frequency (Hz)",ylabel="Amplitude (m/s)",xlog=True,ylog=True,pfa=pfa)
        fig_save(pfa[0],args,ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"spe2")

        #map
        ch_pos_l = np.asarray(ch_pos_list)
        ch_pos_l2 = np.asarray(ch_pos_list_2)
        fig = modules.hypo_mono(ch_pos_l[:,0],ch_pos_l[:,1],c="gray",style="i0.2c",pen="black",lat_range=map_la,lon_range=map_lo)
        fig = modules.hypo_mono(ch_pos_l2[:,0],ch_pos_l2[:,1],c="black",style="i0.2c",pen="black",fig=fig)
        fig = modules.hypo_mono(ch_pos_l2[:,0][s_c_n],ch_pos_l2[:,1][s_c_n],c="orangered",style="i0.2c",pen="black",fig=fig)
        fig.plot(x=lon,y=lat,fill="red",style="a0.3c",pen="black")
        if args.text == "Y":
            text_list = [[map_la[1]-0.1,map_lo[0]+0.2,"%s (s) - %s - %s (s)" %(dur_pre,ymdhms,dur_pos)],
                         [map_la[0]+1.28,map_lo[1]-3,"depth %s (km)" %dep],
                         [map_la[0]+1.08,map_lo[1]-3,"envelope dur %s (s)" %dur],
                         [map_la[0]+0.88,map_lo[1]-3,"energy dur %s (s)" %energy_dur],
                         [map_la[0]+0.68,map_lo[1]-3,"magnitude %s" %df["mag"][i]],
                         [map_la[0]+0.48,map_lo[1]-3,"1st azimuthal gap %s" %df["first az"][i]],
                         [map_la[0]+0.28,map_lo[1]-3,"2nd azimuthal gap %s" %df["second az"][i]],
                         [map_la[0]+0.08,map_lo[1]-3,"ratio %s" %f_ratio]]
            fig = modules.map_text(text_list,font="8p",angle=0,justify="ML",fig=fig)
        fig_save(fig,args,ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"map")

        #scaling
        if args.catalog == "all":
            df_back = pd.read_csv(args.savepath + "/" + args.foldername + "/events.csv")
        else:
            df_back = df
        pfa = modules.simple_plot(df_back["Es"],df_back["energy dur"],
                                  xlog=True,ylog=True,c="darkgray",ec=None)
        pfa[2].set_xlim(10**1,10**10)
        pfa[2].set_ylim(10**-2,10**3)
        pfa = modules.simple_plot(df["Es"][i],energy_dur,xlabel="Seismic energy (J)",ylabel="Minimum energy duration (s)"
                                  ,c="red",s=170,marker="*",ec="black",lw=0.5,pfa=pfa)
        fig_save(pfa[0],args,ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"sca")

        #AH-AL
        pfa = modules.simple_plot(df_back["f low sn"],df_back["f high sn"],xlog=True,ylog=True,c="darkgray",ec=None)
        pfa = modules.simple_plot(f_low_sn,f_high_sn,c="red",xlabel="Mean amplitude at low frequency (m/s)",ylabel="Mean amplitude at high frequency (m/s)",
                                  s=250,marker="*",ec="black",lw=0.5,pfa=pfa)
        pfa[2].set_xlim(10**-9,2*10**-4)
        pfa[2].set_ylim(10**-9,2*10**-4)
        fig_save(pfa[0],args,ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"fre")

        #結合
        if args.outfile == "png":
            map_path = args.savepath + "/eventfig/" + ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"map.png"
            waveform_path = args.savepath + "/eventfig/" + ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"waveform.png"
            pic = modules.get_concat_h_resize(map_path, waveform_path)
            fig_save(pic,args,ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"tent1")
            tent1_path = args.savepath + "/eventfig/" + ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"tent1.png"
            sca_path = args.savepath + "/eventfig/" + ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"sca.png"
            pic = modules.get_concat_h_resize(tent1_path, sca_path)
            fig_save(pic,args,ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"tent2")
            
            spe_path = args.savepath + "/eventfig/" + ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"spe.png"
            spe2_path = args.savepath + "/eventfig/" + ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"spe2.png"
            pic = modules.get_concat_h_resize(spe_path, spe2_path)
            fig_save(pic,args,ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"tent3")
            tent3_path = args.savepath + "/eventfig/" + ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"tent3.png"
            cum_path = args.savepath + "/eventfig/" + ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"cum.png"
            pic = modules.get_concat_h_resize(tent3_path, cum_path)
            fig_save(pic,args,ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"tent4")
            tent4_path = args.savepath + "/eventfig/" + ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"tent4.png"
            fre_path = args.savepath + "/eventfig/" + ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"fre.png"
            pic = modules.get_concat_h_resize(tent4_path, fre_path)
            fig_save(pic,args,ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"tent5")
            
            tent2_path = args.savepath + "/eventfig/" + ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"tent2.png"
            tent5_path = args.savepath + "/eventfig/" + ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"tent5.png"
            pic = modules.get_concat_v_resize(tent2_path, tent5_path)
            fig_save(pic,args,ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i))
            
            for _ in ["map","waveform","sca","spe","spe2","cum","fre","tent1","tent2","tent3","tent4","tent5"]:
                subprocess.call("rm %s" %(args.savepath + "/eventfig/" + ymd+"_"+str(args.ymdhms)[11:19].replace(":","")+"_"+str(i)+"%s.png" %_),shell=True)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--savepath', help='save path',required=True)
    parser.add_argument('--foldername', help='folder',required=True)
    parser.add_argument('--outfile', help='eps output',required=True)
    parser.add_argument('--text', help='text in hypo image',required=True)
    parser.add_argument('--catalog', help='tremors or all events',required=True)
    parser.add_argument('--ymdhms', help='YYYY-MM-DD HH:MM:SS.SSSSSS',required=True)
    args = parser.parse_args()

    print(args.ymdhms)

    if args.catalog == "all":
        df = pd.read_csv(args.savepath + "/" + args.foldername + "/events.csv")
    elif args.catalog == "tremor":
        df = pd.read_csv(args.savepath + "/" + args.foldername + "/tremor_catalog.csv")
    event_fig(args,df) 