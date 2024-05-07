import argparse
import glob
from tqdm import tqdm
import datetime
import pandas as pd
import modules
import numpy as np
import matplotlib.pyplot as plt
import const
import sys


def json_date_check(args,json_path_list):

    all_day = []
    ymd = datetime.datetime(year=int(args.year),month=1,day=1)
    while True:
        all_day.append(args.savepath + "/" + str(ymd.year) + "/json/" + str(ymd.year)[2:4] + str(ymd.month).zfill(2) + str(ymd.day).zfill(2) + ".json")
        ymd += datetime.timedelta(days=1)
        
        if str(ymd.year) != str(args.year):
            break

    lack = False
    for js in all_day:
        if js in json_path_list:
            None
        else:
            print(js,"is not found")
            lack = True
    
    return lack


def originals(args):
    
    print(args.year)
    json_path = args.savepath + "/" + args.year + "/json/" + "%s????.json" %str(args.year)[2:4]
    json_path_list = sorted(glob.glob(json_path))

    lack = json_date_check(args,json_path_list)
    if lack:
        yn = str(input("continue? [y/n] : "))
        if yn == "y":
            None
        elif yn != "y":
            sys.exit()

    data_list = []
    print("json files to csv")
    for j_p in tqdm(json_path_list):
        date = "20" + (j_p.replace(".json","")).replace(args.savepath + "/" + args.year + "/json/","")
        y,m,d = int(date[0:4]), int(date[4:6]), int(date[6:8])
        delta_day = (datetime.datetime(year=y,month=m,day=d) - datetime.datetime(year=int(args.year),month=1,day=1)).days
        df = pd.read_json(j_p)
        
        for k in range(len(df)):
            first_az = float(df.at[k,"1st az"])
            second_az = float(df.at[k,"2nd az"])
            channel_path = df.at[k,"channel path"]
            channel_position = df.at[k,"channel position"]
            cc = float(df.at[k,"cross correlation"])
            dep = float(df.at[k,"depth"])
            dep_error = float(df.at[k,"depth error"])
            dur = int(df.at[k,"duration"])
            dur_pos = int(df.at[k,"duration pos"])
            dur_pre = int(df.at[k,"duration pre"])
            lat = float(df.at[k,"latitude"])
            lat_error = float(df.at[k,"latitude error"])
            lon = float(df.at[k,"longitude"])
            lon_error = float(df.at[k,"longitude error"])
            mag = float(df.at[k,"me"])
            Es = 10**(1.5*mag+4.4)
            doy = df.at[k,"time"]/(3600*24) + delta_day
            ymdhms = datetime.datetime(int(args.year),1,1,0,0,0) + datetime.timedelta(days=doy)

            year, month, day = int(ymdhms.year), int(ymdhms.month), int(ymdhms.day)
            hour, minu, sec, microsec = int(ymdhms.hour), int(ymdhms.minute), int(ymdhms.second), int(ymdhms.microsecond)
         
            a = [ymdhms,year,month,day,hour,minu,sec,microsec,doy,
                 lat,lon,mag,Es,dep,dur,dur_pre,dur_pos,
                 channel_path,channel_position,first_az,second_az,
                 lat_error,lon_error,dep_error,cc]

            data_list.append(a)

    print("saving csv ...")
    df = pd.DataFrame(data_list,columns=["ymdhms","year","month","day","hour","min","sec","microsec","doy",
                                        "lat","lon","mag","Es","dep","dur","dur pre","dur pos",
                                        "channel path","channel pos","first az","second az",
                                        "lat error","lon error","dep error","cc"])
    df.sort_values(by="ymdhms", ascending=True, inplace=True)
    df.reset_index(drop=True, inplace=True)
    df.to_csv(args.savepath + "/" + args.year + "/csv_files/original.csv",index=False)
    print("json to csv : Done")


def removal(df,thres,element_min=2):
    
    c = 1
    cluster_num_list = [c]
    for i in range(1,len(df)):
        if abs(df["scaled dt"][i])<thres[0] and df["scaled dx"][i]<thres[1]:
            cluster_num_list.append(cluster_num_list[int(df["nearest num"][i])])
        else:
            c += 1
            cluster_num_list.append(c)
    cluster_num = []
    for i in range(len(df)):
        if cluster_num_list.count(cluster_num_list[i]) >= element_min:
            cluster_num.append(cluster_num_list[i])
        else:
            cluster_num.append(-1)
    setlist = sorted(list(set(cluster_num)))
    cluster_num = np.array(cluster_num)
    n = 1
    for i in setlist:
        if i == -1:
            None
        else:
            cluster_num[cluster_num==i] = n
            n += 1

    df["clu num"] = cluster_num

    ch_num_list = []
    for i in range(len(df)):
        chs = eval(df["channel path"][i])
        ch_num_list.append(len(chs))

    nums = []
    drop_list = []
    for i in range(len(df)):
        if df["clu num"][i] == -1:
            None
        elif i not in nums:
            df_clu = df[df["clu num"]==df["clu num"][i]].copy()
            clu_index_list = list(df_clu.index)
            clu_ch_num_list = []
            for j in clu_index_list:
                clu_ch_num_list.append(ch_num_list[j])
            max_index = clu_ch_num_list.index(max(clu_ch_num_list))
            for j in range(len(clu_ch_num_list)):
                nums.append(clu_index_list[j])
                if j == max_index:
                    None
                else:
                    drop_list.append(clu_index_list[j])

    df = df.drop(df.index[drop_list])
    df.reset_index(drop=True, inplace=True)    
    return df


def duplication(args):
    
    print("starting duplication removal")
    
    df = pd.read_csv(args.savepath + "/" + args.year + "/csv_files/original.csv")
    print("all detected events in %s (duplication included) :" %args.year,len(df))

    df = modules.nearest(df,"AO",[1,0])
    tf = (0<df["scaled dx"])&(df["scaled dx"]<np.inf)&(0<df["scaled dt"])&(df["scaled dt"]<np.inf)
    
    hist, bin_edges = np.histogram(df[tf]["scaled dt"],bins=np.logspace(-9,1,50))
    bin_m = []
    for b in range(len(bin_edges)-1):
        bin_m.append((bin_edges[b]+bin_edges[b+1])/2)
    bin_m = np.array(bin_m)
    peak1,peak2 = max(hist[bin_m<10/86400]),max(hist[bin_m>10/86400])
    bin1,bin2 = bin_m[hist==peak1],bin_m[hist==peak2]
    _,bin = max(hist),None

    for j in range(len(hist)):
        if bin1[0]<bin_m[j] and bin_m[j]<bin2[0] and hist[j]<_:
            _,bin = hist[j],bin_m[j]

    pfa = modules.simple_plot(df[tf]["scaled dt"],df[tf]["scaled dx"],xlabel="ΔT (day)",ylabel="ΔR (km)",xlog=True,ylog=True,s=1,marker=".")
    t = np.array([10**-20,10**20])
    pfa[2].plot(t,t*const.beta/10**3*86400,linestyle="dashed",c="black",label="Vs=%s (km/s)" %str(const.beta/10**3))
    pfa[2].set_xlim(10**-9, 10**1)
    pfa[2].set_ylim(10**-3, 10**3)
    pfa[2].plot([bin,bin],[10**-6,10**4],linewidth=0.5,c="black",linestyle="dashed")
    pfa[2].plot([10**-10,10**2],[bin*const.beta/10**3*86400]*2,linewidth=0.5,c="black",linestyle="dashed")
    pfa[2].legend()
    pfa[0].savefig(args.savepath + "/" + args.year + "/figs/duplication.png")

    print("time difference of duplication is automatically determined as : ",bin*86400,"(sec)")
    _ = str(input("Do you use this value as threshold? [y/n] : "))
    print("removing duplications ...")
    if _ == "y":
        th = [bin, bin*const.beta/10**3*86400]
    else:
        delta_t = float(input("input time difference (day): "))
        th = [delta_t, delta_t*const.beta/10**3*86400]

    df = removal(df,th,element_min=2)

    df.drop(['scaled dx','scaled dt','eta','nearest num','nearest ymdhms','nearest mag','nearest dep','clu num'], axis=1, inplace=True)
    df = df[df["first az"]<const.azimuthal_gap].copy()
    df.sort_values(by="ymdhms", ascending=True, inplace=True)
    df.reset_index(drop=True, inplace=True)
    df.to_csv(args.savepath + "/" + args.year + "/csv_files/events.csv",index=False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--year', help='year',required=True)
    parser.add_argument('--savepath', help='save path',required=True)
    args = parser.parse_args()
    
    originals(args)
    duplication(args)