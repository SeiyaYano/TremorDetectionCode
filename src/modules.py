import pyproj
import numpy as np
from scipy.interpolate import RectBivariateSpline
import pygmt
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy import interpolate
import matplotlib as mpl
plt.rcParams["font.size"] = 13
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True

def inter(x,y,kind,x0):
    f = interpolate.interp1d(x, y, kind="%s" %kind, bounds_error=True)
    return f(x0)


def mean_norm(df_input,axis=0):
    return df_input.apply(lambda x: (x-x.mean())/ x.std(), axis=axis)


def latlon2xy(origin,point):
    #originを中心としたpointの方位角[degree](北が0度, 時計回りが正)と距離[km]
    grs80 = pyproj.Geod(ellps='GRS80')
    azi, bkw_azi, dis = grs80.inv(origin[1],origin[0],point[1],point[0])
    dis = dis/1000.0    #[km]
    x, y = dis*np.sin(azi*np.pi/180), dis*np.cos(azi*np.pi/180)
    return {
        "azimuth" : azi,
        "distance" : dis,
        "x" : x,
        "y" : y
    }

def japan_land_travel_time(hypo_pos,sta_pos,P_or_S):
    #hypo_pos; [lat,lon,dep]
    if len(hypo_pos) != 3:
        print("confirm that hypo_pos has lat, lon, and dep")
    if P_or_S == "P":
        x, y = [], []
        t = np.empty((106, 236), dtype=np.float32)
        with open("JMA2001/travel_time") as f:
            for line in f.readlines():
                _, pt, _, _, dep, dis = line.split()
                pt, dep, dis = map(float, (pt, dep, dis))
                if dep not in x:
                    x.append(dep)
                if dis not in y:
                    y.append(dis)
                t[x.index(dep), y.index(dis)] = pt
        _time = RectBivariateSpline(x, y, t, kx=1, ky=1)
        dis = latlon2xy(sta_pos,hypo_pos)["distance"]
        dt = _time.ev(hypo_pos[2],dis)
    if P_or_S == "S":
        x, y = [], []
        t = np.empty((106, 236), dtype=np.float32)
        with open("JMA2001/travel_time") as f:
            for line in f.readlines():
                _, _, _, st, dep, dis = line.split()
                st, dep, dis = map(float, (st, dep, dis))
                if dep not in x:
                    x.append(dep)
                if dis not in y:
                    y.append(dis)
                t[x.index(dep), y.index(dis)] = st
        _time = RectBivariateSpline(x, y, t, kx=1, ky=1)
        dis = latlon2xy(sta_pos,hypo_pos)["distance"]
        dt = _time.ev(hypo_pos[2],dis)
    return float(dt)

def hypo_mono(lat,lon,c="orangered",style="c0.07c",pen="black",sl="0.1p,gray",lat_range=None,lon_range=None,fig=None):
    if fig is None:
        fig = pygmt.Figure()
        if lat_range is None and lon_range is None:
            fig.basemap(region=[min(lon)-0.25,max(lon)+0.25,min(lat)-0.25,max(lat)+0.25],projection="M12c",frame="a")
        else:
            fig.basemap(region=[lon_range[0],lon_range[1],lat_range[0],lat_range[1]],projection="M12c",frame="a")
        fig.coast(land="lightgray",water="white",shorelines=sl)
    fig.plot(x=lon,y=lat,fill=c,style=style,pen=pen)
    return fig

def hypo_colorbar(lat,lon,c,clabel,clog=False,cpos=None,style="c0.07c",pen="0.05p,black",sh="0.1p,gray",
                  lat_range=None,lon_range=None,crange=None,fig=None):
    if fig is None:
        fig = pygmt.Figure()
        if lat_range is None and lon_range is None:
            delta_lat = (max(lat)-min(lat))*0.05
            delta_lon = (max(lon)-min(lon))*0.05
            fig.basemap(region=[min(lon)-delta_lon,max(lon)+delta_lon,min(lat)-delta_lat,max(lat)+delta_lat],projection="M12c",frame="a")
        else:
            fig.basemap(region=[lon_range[0],lon_range[1],lat_range[0],lat_range[1]],projection="M12c",frame="a")
        fig.coast(land="lightgray",water="white",shorelines=sh)
    if crange is None:
        if not clog:
            delta_c = (max(c)-min(c))*0.05
            pygmt.makecpt(cmap="jet",series=[min(c)-delta_c, max(c)+delta_c],log=clog,no_bg=clog)
        else:
            pygmt.makecpt(cmap="jet",series=[np.log10(min(c))-0.2,np.log10(max(c))+0.2],log=clog,no_bg=clog)
    else:
        if not clog:
            delta_c = (crange[1]-crange[0])*0.05
            pygmt.makecpt(cmap="jet",series=[crange[0]-delta_c, crange[1]+delta_c],log=clog,no_bg=clog)
        else:
            pygmt.makecpt(cmap="jet",series=[np.log10(crange[0])-0.2,np.log10(crange[1])+0.2],log=clog,no_bg=clog)
    fig.plot(x=lon,y=lat,fill=c,cmap=True,style=style,pen=pen)
    fig.colorbar(frame="af+l%s"%clabel,log=clog,position=cpos)
    return fig

def plate_contour(pathlist,pen="0.5,black,-",sl="0.1p,gray",lat_range=None,lon_range=None,fig=None):
    if fig is None:
        fig = pygmt.Figure()
        fig.basemap(region=[lon_range[0],lon_range[1],lat_range[0],lat_range[1]],projection="M12c",frame="a")
        fig.coast(land="lightgray",water="white",shorelines=sl)
    for path in pathlist:
        contour = pd.read_csv(path)
        fig.plot(x=contour["lon"],y=contour["lat"],pen=pen)
    return fig

def volcano(df_vol,c="red",style="t0.2c",pen="0.1p,black",tp=20,sl="0.1p,gray",lat_range=None,lon_range=None,fig=None):
    if fig is None:
        fig = pygmt.Figure()
        fig.basemap(region=[lon_range[0],lon_range[1],lat_range[0],lat_range[1]],projection="M12c",frame="a")
        fig.coast(land="lightgray",water="white",shorelines=sl)
    fig.plot(x=df_vol["lon"],y=df_vol["lat"],fill=c,style=style,pen=pen,transparency=tp)
    return fig

def map_text(pos_text_set_list,font="10p",angle=0,sl="0.1p,gray",lat_range=None,lon_range=None,justify=None,fig=None):
    #pos_text_set_list = [[lat,lon,text],[lat,lon,text],⋯]
    if fig is None:
        fig = pygmt.Figure()
        fig.basemap(region=[lon_range[0],lon_range[1],lat_range[0],lat_range[1]],projection="M12c",frame="a")
        fig.coast(land="lightgray",water="white",shorelines=sl)
    for pos_text_set in pos_text_set_list:
        fig.text(x=pos_text_set[1],y=pos_text_set[0],text=pos_text_set[2],font=font,angle=angle,justify=justify)
    return fig

def get_concat_h_resize(im1_path, im2_path, resample=Image.BICUBIC, resize_big_image=True):
    im1 = Image.open(im1_path)
    im2 = Image.open(im2_path)
    if im1.height == im2.height:
        _im1 = im1
        _im2 = im2
    elif (((im1.height > im2.height) and resize_big_image) or
          ((im1.height < im2.height) and not resize_big_image)):
        _im1 = im1.resize((int(im1.width * im2.height / im1.height), im2.height), resample=resample)
        _im2 = im2
    else:
        _im1 = im1
        _im2 = im2.resize((int(im2.width * im1.height / im2.height), im1.height), resample=resample)
    dst = Image.new('RGB', (_im1.width + _im2.width, _im1.height))
    dst.paste(_im1, (0, 0))
    dst.paste(_im2, (_im1.width, 0))
    return dst

def get_concat_v_resize(im1_path, im2_path, resample=Image.BICUBIC, resize_big_image=True):
    im1 = Image.open(im1_path)
    im2 = Image.open(im2_path)
    if im1.width == im2.width:
        _im1 = im1
        _im2 = im2
    elif (((im1.width > im2.width) and resize_big_image) or
          ((im1.width < im2.width) and not resize_big_image)):
        _im1 = im1.resize((im2.width, int(im1.height * im2.width / im1.width)), resample=resample)
        _im2 = im2
    else:
        _im1 = im1
        _im2 = im2.resize((im1.width, int(im2.height * im1.width / im2.width)), resample=resample)
    dst = Image.new('RGB', (_im1.width, _im1.height + _im2.height))
    dst.paste(_im1, (0, 0))
    dst.paste(_im2, (0, _im1.height))
    return dst


def nearest(df,eta_def,search_range,dur_def=None,dur_unit=None,f_d=None):

    if dur_def is None or dur_def not in df.columns:
        dur_def = "tentative"
        df[dur_def] = -1
    lat,lon = df["lat"],df["lon"]
    mag,dur,dep,ymdhms = df["mag"],df[dur_def],df["dep"],df["ymdhms"]
    doy = df["doy"]

    if dur_unit == "m":
        dur = dur * 60
    elif dur_unit == "h":
        dur = dur * 60*60
    elif dur_unit == "d":
        dur = dur * 60*60*24
    dx_list,dt_list,eta_list = [],[],[]
    p_num,p_ymdhms,p_dur,p_mag,p_dep = [],[],[],[],[]
    group = [dx_list,dt_list,eta_list,p_num,p_ymdhms,p_dur,p_mag,p_dep]
    
    for i in tqdm(range(len(df))):
        df_i = df[(df["doy"][i]-search_range[0]<df["doy"]) & (df["doy"]<df["doy"][i]+search_range[1])].copy()
        info = [np.inf,np.inf,np.inf,None,None,None,None,None]
        for j in df_i.index.values:
            if j == i:
                continue
            dis = latlon2xy([lat[i],lon[i]],[lat[j],lon[j]])["distance"]
            dis = ((dis**2 + (dep[i]-dep[j])**2)**0.5)  #[km]
            dt = (doy[i] - doy[j])  #[day]
            if eta_def == "AO":
                dx,dt = dis,dt
            elif eta_def == "ZB":
                dis = dis ** f_d
                dx,dt = dis*10**(-0.5*mag[j]),dt*10**(-0.5*mag[j])
            elif eta_def == "mZB":
                dis = dis ** f_d
                dx, dt = dis*10**(-0.5*mag[j]),  dt/dur[j]
            if abs(dx*dt) < info[2]:
                info[0],info[1],info[2] = dx,dt,dx*dt
                info[3],info[4],info[5] = j,ymdhms[j],dur[j]
                info[6],info[7] = mag[j],dep[j]
        for p in range(len(group)):
            group[p].append(info[p])
    df["scaled dx"],df["scaled dt"],df["eta"] = dx_list, dt_list, eta_list
    df["nearest num"],df["nearest ymdhms"] = p_num, p_ymdhms
    df["nearest dur"],df["nearest mag"],df["nearest dep"] = p_dur, p_mag, p_dep
    if df[dur_def][0] < 0:
        df.drop(columns="tentative", inplace=True)
        df.drop(columns='nearest dur', inplace=True)
    return df

def simple_plot(x,y,legend=None,xlabel=None,ylabel=None,xlog=False,ylog=False,c="black",s=10,marker=".",
                ec="black",alpha=1,lw=0.1,figsize=None,pfa=[None]*3):

    if pfa is None:
        pfa = [None]*3

    if pfa[0] is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
    elif pfa[0] is not None:
        fig = pfa[1]
        ax = pfa[2]
    
    ax.scatter(x,y,c=c,s=s,marker=marker,ec=ec,lw=lw,alpha=alpha,label=legend)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    if xlog:
        ax.set_xscale('log')
    if ylog:
        ax.set_yscale('log')
    
    if legend is not None:
        ax.legend()
    
    return (plt,fig,ax)


def simple_graph(x,y,xlabel=None,ylabel=None,legend=None,xlog=False,ylog=False,c="black",alpha=1,lw=1,figsize=None,pfa=[None]*3):

    if pfa is None:
        pfa = [None]*3

    if pfa[0] is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
    elif pfa[0] is not None:
        fig = pfa[1]
        ax = pfa[2]
    
    ax.plot(x,y,c=c,lw=lw,alpha=alpha,label=legend)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    if xlog:
        ax.set_xscale('log')
    if ylog:
        ax.set_yscale('log')
    
    if legend is not None:
        ax.legend()
    
    return (plt,fig,ax)


def simple_histogram(x,bins,label,legend=None,xrange=None,xlog=False,ylog=False,cum=False,
                     c="cornflowerblue",alpha=1.0,ec="black",figsize=None,pfa=[None]*3):
    
    if pfa is None:
        pfa = [None]*3

    if pfa[0] is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
    elif pfa[0] is not None:
        fig = pfa[1]
        ax = pfa[2]
    
    if xrange is None:
        if xlog:
            bins = np.logspace(np.log10(min(x)-min(x)*0.05),np.log10(max(x)+max(x)*0.05),bins)
            ax.set_xscale('log')
        elif not xlog:
            delta = (max(x)-min(x))*0.05
            bins = np.linspace(min(x)-delta,max(x)+delta,bins)
    elif xrange is not None:
        if xlog:
            bins = np.logspace(np.log10(xrange[0])-0.2,np.log10(xrange[1])+0.2,bins)
            ax.set_xscale("log")
        elif not xlog:
            delta = (xrange[1] - xrange[0])*0.05
            bins = np.linspace(xrange[0]-delta,xrange[1]+delta,bins)
    
    ax.hist(x,bins=bins,color=c,histtype="bar",alpha=alpha,ec=ec,log=ylog,cumulative=cum,label=legend)
    ax.set_xlabel(label)
    ax.set_ylabel("Number")

    if legend is not None:
        ax.legend()
    
    return (plt,fig,ax)


def two_d_histogram(x,xlabel,y,ylabel,bins=100,log=False,clog=False,clabel=None):
    if log:
        x = x[(x>0)&(y>0)].copy()
        y = y[(x>0)&(y>0)].copy()
    maxi, mini = max(max(x),max(y)), min(min(x),min(y))
    fig, ax = plt.subplots()
    if log and clog:
        H = ax.hist2d(x,y,bins=np.logspace(np.log10(mini)-0.2,np.log10(maxi)+0.2,bins),cmap="Blues",weights=np.ones_like(x),norm=mpl.colors.LogNorm())
    elif not log and clog:
        delta = (maxi - mini)*0.05
        H = ax.hist2d(x,y,bins=np.linspace(mini-delta,maxi+delta,bins),cmap="Blues",weights=np.ones_like(x),norm=mpl.colors.LogNorm())
    elif log and not clog:
        H = ax.hist2d(x,y,bins=np.logspace(np.log10(mini)-0.2,np.log10(maxi)+0.2,bins),cmap="Blues",weights=np.ones_like(x))
    elif not log and not clog:
        delta = (maxi - mini)*0.05
        H = ax.hist2d(x,y,bins=np.linspace(mini-delta,maxi+delta,bins),cmap="Blues",weights=np.ones_like(x))
    HH = H[0]
    if clog:
        HH[HH == 0] = 1
        HH = np.log10(HH)
    grid = HH.transpose()
    grid_padded = np.pad(grid,[(0,1)])
    ax.contourf(H[1],H[2],grid_padded,cmap="Blues")
    if log:
        ax.set_xscale('log')
        ax.xaxis.set_major_locator(mpl.ticker.LogLocator(numticks=50))
        ax.set_yscale('log')
        ax.yaxis.set_major_locator(mpl.ticker.LogLocator(numticks=50))
        ax.xaxis.set_minor_locator(mpl.ticker.LogLocator(numticks=50, subs=(.1, .2, .3, .4, .5, .6, .7, .8, .9)))
        ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(numticks=50, subs=(.1, .2, .3, .4, .5, .6, .7, .8, .9)))
    ax.set_xlabel("%s" % xlabel)
    ax.set_ylabel("%s" % ylabel)
    ax.grid(which="major", axis="x", color="gray", alpha=0.3, linestyle="--", linewidth=0.6)
    ax.grid(which="major", axis="y", color="gray", alpha=0.3, linestyle="--", linewidth=0.6)
    ax.grid(which="minor", axis="x", color="gray", alpha=0.3, linestyle="--", linewidth=0.6)
    ax.grid(which="minor", axis="y", color="gray", alpha=0.3, linestyle="--", linewidth=0.6)
    cbar = fig.colorbar(H[3],ax=ax)
    cbar.set_label(clabel)
    return (plt,fig,ax)

