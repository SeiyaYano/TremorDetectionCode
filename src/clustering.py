import argparse
import algorithm
import glob
import pandas as pd
import modules as md
import const
import numpy as np
import subprocess
from scipy.stats import boxcox
import sys
from tqdm import tqdm
from decimal import *
import matplotlib.pyplot as plt
import os
from concurrent.futures import ProcessPoolExecutor
import dimension
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
import datetime
import GPy, GPyOpt


def events(args):
    
    print("preparing ...")
    
    if args.year == "all":
        csv_path = args.savepath + "/????/csv_files/features.csv"
    else:
        csv_path = args.savepath + "/" + args.year + "/csv_files/features.csv"
    csv_list = sorted(glob.glob(csv_path))

    dt_now = datetime.datetime.now()
    f = open(args.savepath + "/" + args.foldername + "/" + "parameters.txt", "a")
    f.write("\n%s\n" %dt_now)
    f.write("%s\n" %str(csv_list))
    f.close()

    for _ in csv_list:
        print(_)

    df_list = []
    for c in csv_list:
        df_p = pd.read_csv(c)
        df_list.append(df_p)

    df = pd.concat(df_list)
    df.sort_values(by ="ymdhms", ascending = True, inplace = True)
    df.reset_index(drop=True, inplace=True)
    df.to_csv(args.savepath + "/" + args.foldername + "/events.csv",index=False)

    f = open(args.savepath + "/" + args.foldername + "/" + "parameters.txt", "a")
    f.write('# of all detected events (azi≦180) : %s\n' %str(len(df)))
    f.close()

    #epicenter
    fig = md.hypo_colorbar(df["lat"],df["lon"],df["dep"],"Depth (km)",clog=True,
                           lat_range=const.lat_range,lon_range=const.lon_range)
    fig.savefig(args.savepath + "/" + args.foldername + "/figs/all_epi.png")
    
    #depth
    pfa = md.simple_histogram(df[df["dep"]<=100]["dep"],50,"Depth (km)")
    pfa[0].savefig(args.savepath + "/" + args.foldername + "/figs/all_dep.png",bbox_inches="tight",pad_inches=0.05)
    pfa[0].clf()
    pfa[0].close()
    
    #scaling
    pfa = md.two_d_histogram(df["Es"],"Seismic energy (J)",
                             df["energy dur"],"Minimum energy duration (s)",log=True)
    pfa[2].set_xlim(10**1,10**6)
    pfa[2].set_ylim(10**-2,10**2)
    pfa[0].savefig(args.savepath + "/" + args.foldername + "/figs/all_sca.png",bbox_inches="tight",pad_inches=0.05)
    pfa[0].clf()
    pfa[0].close()

    #AH-AL
    pfa = md.two_d_histogram(df["f low sn"],"Mean amplitude at low frequency (m/s)",
                             df["f high sn"],"Mean amplitude at high frequency (m/s)",log=True)
    pfa[2].set_xlim(10**-8.5,10**-5)
    pfa[2].set_ylim(10**-8.5,10**-5)
    x = np.array([10**-10,10**10])
    pfa[2].plot(x,x,label="AH/AL = 1.0",linestyle="dashed",color="black",linewidth=0.7)
    pfa[2].plot(x,x*0.4,label="AH/AL = 0.4",linestyle="dashed",color="gray",linewidth=0.7)
    pfa[2].plot(x,x*0.2,label="AH/AL = 0.2",linestyle="dashed",color="lightgray",linewidth=0.7)
    pfa[2].legend()
    pfa[0].savefig(args.savepath + "/" + args.foldername + "/figs/all_fre.png",bbox_inches="tight",pad_inches=0.05)
    pfa[0].clf()
    pfa[0].close()

    #dur
    pfa = md.simple_plot(df["dur"],df["energy dur"],
                         xlabel="Envelope duration (s)",ylabel="Minimum energy duration (s)",
                         xlog=True,ylog=True)
    pfa[2].set_xlim(1,10**3)
    pfa[2].set_ylim(0.01,10**3)
    x = np.array([0,10**4])
    pfa[2].plot(x,x,lw=1,c="black")
    pfa[2].plot(x,x/2,lw=1,c="black",linestyle="dashed")
    pfa[0].savefig(args.savepath + "/" + args.foldername + "/figs/all_dur.png",bbox_inches="tight",pad_inches=0.05)
    pfa[0].clf()
    pfa[0].close()


class Clustering():
    
    def __init__(self,args):
        self.args = args
        columns_list = ["ymdhms","year","month","day","hour","min","sec","microsec","doy",
                        "lat","lon","mag","Es","dep",
                        "dur","dur pre","dur pos",
                        "channel path","channel pos",
                        "first az","second az",
                        "lat error","lon error","dep error",
                        "cc",
                        "energy dur","energy dur list",
                        "f low e","f high e","f low n","f high n","f low sn","f high sn",
                        "f ratio",
                        "channel pos 2"]
        self.columns = columns_list
        

    def clustering(self):

        df = pd.read_csv(self.args.savepath + "/" + self.args.foldername + "/events.csv")
        df_features = pd.DataFrame()
        
        df_features["energy dur"] = boxcox(df["energy dur"])[0]
        df_features["dep"] = boxcox(df["dep"])[0]
        df_features["f ratio"] = boxcox(df["f ratio"])[0]
        df_features["f high sn"] = boxcox(df["f high sn"])[0]
        df_features["f low sn"] = boxcox(df["f low sn"])[0]

        if const.min_s == 0:
            print("searching for the best p")
            if len(df) >= 10000:
                p_range = [0.05,1.0]
            elif len(df) >= 1000:
                p_range = [0.5,10.0]
            elif len(df) < 1000:
                print("require more than 1000 detected events")
                sys.exit()

            min_s_list = []
            score_list = []
            # for min_s in tqdm(range(int(len(df)*p_range[0]),int(len(df)*p_range[1])+1)):
            #     score = algorithm.HDBSCAN(df_features,min_s)["score"]
            #     min_s_list.append(min_s)
            #     score_list.append(score)
            
            u,l = int(len(df)*p_range[1]/100)+1,int(len(df)*p_range[0]/100)
            with tqdm(total=u-l) as progress:
                with ProcessPoolExecutor(max_workers=3) as executor:
                    for min_s in range(l,u):
                        future = executor.submit(algorithm.HDBSCAN,df_features,min_s)
                        future.add_done_callback(lambda p: progress.update())
                        min_s_list.append(min_s)
                        res = future.result()
                        score_list.append(res["score"])

            plt.figure()
            plt.plot(min_s_list,score_list,lw=1.5)
            plt.xscale("log")
            plt.xlabel("minimum cluster size")
            plt.ylabel("score")
            plt.savefig(self.args.savepath + "/" + self.args.foldername + "/clustering/p_values.png")
            plt.clf()
            plt.close()

            min_s_best = min_s_list[score_list.index(max(score_list))]
        else:
            min_s_best = const.min_s
            
        res = algorithm.HDBSCAN(df_features,min_s_best)

        print("min # of events : ",min_s_best," / # of clusters : ",max(res["clu num"])+1)
        f = open(args.savepath + "/" + args.foldername + "/" + "parameters.txt", "a")
        f.write('minimum # of events: %s\n' %str(min_s_best))
        f.write('# of clusters: %s\n' %str(max(res["clu num"])+1))
        f.close()

        df["clu num"] = res["clu num"]
        
        c_num = sorted(list(set(df["clu num"])))
        for i in c_num:
            df_i = df[df["clu num"]==i].copy()
            df_i.sort_values(by ="ymdhms", ascending = True, inplace = True)
            df_i.reset_index(drop=True, inplace=True)
            
            #epicenter
            fig = md.hypo_colorbar(df_i["lat"],df_i["lon"],df_i["dep"],"Depth (km)",
                                lat_range=const.lat_range,lon_range=const.lon_range)
            fig.savefig(self.args.savepath + "/" + self.args.foldername + "/clustering/c%s_epi.png" %i)
            
            #depth
            pfa = md.simple_histogram(df_i["dep"],50,"Depth (km)")
            pfa[0].savefig(self.args.savepath + "/" + self.args.foldername + "/clustering/c%s_dep.png" %i,bbox_inches="tight",pad_inches=0.05)
            pfa[0].clf()
            pfa[0].close()

            #AH-AL
            pfa = md.simple_plot(df_i["f low sn"],df_i["f high sn"],xlabel="AL (m/s)",ylabel="AH (m/s)",
                                xlog=True,ylog=True)
            pfa[2].set_xlim(10**-8.5,10**-5)
            pfa[2].set_ylim(10**-8.5,10**-5)
            x = np.array([10**-10,10**10])
            pfa[2].plot(x,x,label="AH/AL = 1.0",linestyle="dashed",color="black",linewidth=0.7)
            pfa[2].plot(x,x*0.4,label="AH/AL = 0.4",linestyle="dashed",color="gray",linewidth=0.7)
            pfa[2].plot(x,x*0.2,label="AH/AL = 0.2",linestyle="dashed",color="lightgray",linewidth=0.7)
            pfa[2].legend()
            pfa[0].savefig(self.args.savepath + "/" + self.args.foldername + "/clustering/c%s_fre.png" %i,bbox_inches="tight",pad_inches=0.05)
            pfa[0].clf()
            pfa[0].close()

            #scaling
            pfa = md.simple_plot(df_i["Es"],df_i["energy dur"],
                                xlabel="Seismic energy (J)",ylabel="Minimum energy duration (s)",
                                xlog=True,ylog=True)
            pfa[2].set_xlim(10**1,10**6)
            pfa[2].set_ylim(10**-2,10**2)
            pfa[0].savefig(self.args.savepath + "/" + self.args.foldername + "/clustering/c%s_sca.png" %i,bbox_inches="tight",pad_inches=0.05)
            pfa[0].clf()
            pfa[0].close()

        tremor_catalogs = []
        name_list,tremors = [],[]
        for clu_num in c_num:
            savename = str(input("label the cluster %s : " %clu_num))
            yn = str(input("tremors ? [y/n] : "))
            name_list.append(savename)
            
            catalog_i = df[df["clu num"]==clu_num].copy()
            catalog_i = catalog_i[self.columns].copy()
            catalog_i.sort_values(by ="ymdhms", ascending = True, inplace = True)
            catalog_i.reset_index(drop=True, inplace=True)
            catalog_i.to_csv(self.args.savepath + "/" + self.args.foldername + "/clustering/c_%s.csv"%savename,index=False)

            if yn == "y":
                tremors.append(savename)
                tremor_catalogs.append(catalog_i)

        self.name_list = name_list
        self.tremors = tremors

        self.tremor_catalog = pd.concat(tremor_catalogs)
        
        if const.fractal_dimenstion == 0:
            dim = dimension.box_counting(self.tremor_catalog)
        else:
            dim = const.fractal_dimenstion
        self.dim = dim
        print("fractal dimension : ",dim)
        f = open(args.savepath + "/" + args.foldername + "/" + "parameters.txt", "a")
        f.write('fractal dimension : %s\n' %str(dim))
        f.close()

    def compile_AB(self):

        if const.eta > 0 and const.thre_T > 0 and const.thre_R > 0:
            tremor_catalog = pd.read_csv(self.args.savepath + "/" + self.args.foldername + "/nearest.csv")
        else:
            tremor_catalog = md.nearest(self.tremor_catalog,"mZB",[3.5,3.5],
                                        dur_def="energy dur",dur_unit="s",f_d=self.dim)
            tremor_catalog.to_csv(self.args.savepath + "/" + self.args.foldername + "/nearest.csv",index=False)

        pfa = md.two_d_histogram(abs(tremor_catalog["scaled dt"]),"scaled T",
                                abs(tremor_catalog["scaled dx"]),"scaled R",log=True,clog=True)
        x = np.logspace(-20,20,50)
        eta_l = [0.01,0.1,1.0,10,100]
        line_color = ["red","coral", "tan","forestgreen","dodgerblue"]
        for i in range(len(eta_l)):
            pfa[2].plot(x,eta_l[i]/x,alpha=0.5,color=line_color[i],label="eta = %s" %eta_l[i])
        pfa[2].legend()
        pfa[2].set_xlim(10**-8, 10**3)
        pfa[2].set_ylim(10**-4, 10**7)
        pfa[0].savefig(self.args.savepath + "/" + self.args.foldername + "/clustering/mZB.png")
        pfa[0].clf()
        pfa[0].close()

        pfa = md.simple_histogram(abs(tremor_catalog["eta"]),50,"eta",xlog=True,ylog=True)
        pfa[0].savefig(self.args.savepath + "/" + self.args.foldername + "/clustering/mZB_eta.png")
        pfa = md.simple_histogram(abs(tremor_catalog["scaled dx"]),50,"scaled R",xlog=True,ylog=True)
        pfa[0].savefig(self.args.savepath + "/" + self.args.foldername + "/clustering/mZB_R.png")
        pfa = md.simple_histogram(abs(tremor_catalog["scaled dt"]),50,"scaled T",xlog=True,ylog=True)
        pfa[0].savefig(self.args.savepath + "/" + self.args.foldername + "/clustering/mZB_T.png")
        pfa[0].clf()
        pfa[0].close()

        def AB():
            if const.eta == 0 or const.thre_T == 0 or const.thre_R == 0:
                ths = []
                _ = float(input("threshold for eta: "))
                ths.append(_)
                _ = float(input("threshold for T: "))
                ths.append(_)
                _ = float(input("threshold for R: "))
                ths.append(_)
            else:
                ths = [const.eta,const.thre_T,const.thre_R]
            
            clustered = tremor_catalog[(tremor_catalog["scaled dt"]<=ths[1]) & (tremor_catalog["scaled dx"]<=ths[2]) & (tremor_catalog["eta"]<=ths[0])].copy()
            isolated = tremor_catalog[(tremor_catalog["scaled dt"]>ths[1]) | (tremor_catalog["scaled dx"]>ths[2]) | (tremor_catalog["eta"]>ths[0]) ].copy()

            clustered_catalog = pd.DataFrame()
            isolated_catalog = pd.DataFrame()
            for c in self.columns:
                clustered_catalog[c] = clustered[c]
                isolated_catalog[c] = isolated[c]

            clustered_catalog.sort_values(by ="ymdhms", ascending = True, inplace = True)
            clustered_catalog.reset_index(drop=True, inplace=True)
            clustered_catalog.to_csv(self.args.savepath + "/" + self.args.foldername + "/clustering/QualityA.csv",index=False)
            
            fig = md.hypo_colorbar(clustered_catalog["lat"],clustered_catalog["lon"],
                                   clustered_catalog["dep"],"Depth (km)",
                                   lat_range=const.lat_range,lon_range=const.lon_range)
            fig.savefig(self.args.savepath + "/" + self.args.foldername + "/clustering/QualityA.png")
            
            isolated_catalog.sort_values(by ="ymdhms", ascending = True, inplace = True)
            isolated_catalog.reset_index(drop=True, inplace=True)
            isolated_catalog.to_csv(self.args.savepath + "/" + self.args.foldername + "/clustering/isolated.csv",index=False)
            
            fig = md.hypo_colorbar(isolated_catalog["lat"],isolated_catalog["lon"],
                                isolated_catalog["dep"],"Depth (km)",
                                lat_range=const.lat_range,lon_range=const.lon_range)
            fig.savefig(self.args.savepath + "/" + self.args.foldername + "/clustering/isolated.png")
        
            return ths
        
        while True:
            ths = AB()
            reset = str(input("finish compiling quality A catalog ? [y/n] : "))
            if reset == "y":
                f = open(args.savepath + "/" + args.foldername + "/" + "parameters.txt", "a")
                f.write('threthod for eta : %s\n' %str(ths[0]))
                f.write('threthod for T : %s\n' %str(ths[1]))
                f.write('threthod for R : %s\n' %str(ths[2]))
                f.close()
                break
            else:
                print("reset thresholds for T, R, and eta")

    def reclassification(self):

#        cols = ["energy dur","f low sn","f high sn","f ratio","dep"]
        cols = ["energy dur","f low sn","f high sn","f ratio"]

        df_all = pd.read_csv(self.args.savepath + "/" + self.args.foldername + "/events.csv")

        f = open(args.savepath + "/" + args.foldername + "/" + "parameters.txt", "a")
        f.write('features for reclassification (neural network) : %s\n' %str(cols))
        f.close()
        
        df_clu_list = []
        i = 1
        for clu in self.name_list:
            if clu != "noise" and clu not in self.tremors:
                df_clu = pd.read_csv(self.args.savepath + "/" + self.args.foldername + "/clustering/c_%s.csv" %clu)
                df_clu["label"] = i
                df_clu_list.append(df_clu)
                i += 1

        tremorA = pd.read_csv(self.args.savepath + "/" + self.args.foldername + "/clustering/QualityA.csv")
        tremor_clu_num = i
        tremorA["label"] = tremor_clu_num
        df_clu_list.append(tremorA)
        
        ans_raw = pd.concat(df_clu_list)
        
        ans_input = pd.DataFrame()
        for col in cols:
            ave_all = np.average(df_all[col])
            std_all = np.average(df_all[col])
            ans_input[col] = (ans_raw[col] - ave_all) / std_all
        
        #noise event
        noise1 = pd.read_csv(self.args.savepath + "/" + self.args.foldername + "/clustering/c_noise.csv")
        noise2 = pd.read_csv(self.args.savepath + "/" + self.args.foldername + "/clustering/isolated.csv")
        noise_raw = pd.concat([noise1,noise2])

        depth_shallow = float(input("minimum depth [km] : "))
        depth_deep = float(input("maximum depth [km] : "))
        noise_raw = noise_raw[(depth_shallow<=noise_raw["dep"]) & (noise_raw["dep"]<=depth_deep)].copy()
        
        f = open(args.savepath + "/" + args.foldername + "/" + "parameters.txt", "a")
        f.write('minimum depth : %s\n' %str(depth_shallow))
        f.write('maximum depth : %s\n' %str(depth_deep))
        f.close()
        
        noise_raw.reset_index(drop=True,inplace=True)
        
        noise_input = pd.DataFrame()
        for col in cols:
            ave_all = np.average(df_all[col])
            std_all = np.average(df_all[col])
            noise_input[col] = (noise_raw[col] - ave_all) / std_all

        #hypoer parameter tuning
        print("searching for the best hyper parameters of NN")
        print("neuron,layer,penalty,activation -> score ave, score std")
        
        parameter_set = []
        for n_ in [3,4,5,6]:
            for l_ in [1,2,3,4,5]:
                parameter_set.append([n_,l_])
        
        n_best,l_best = 0,0
        s_best = 0
        for ps in parameter_set:
            neuron = ps[0]
            layer = ps[1]
            
            score_list = []
            for r in [0,1,2,3,4]:
                data_train,data_test,target_train,target_test = train_test_split(ans_input,ans_raw["label"],
                                                                                 test_size=0.4,random_state=r,stratify=ans_raw["label"])
                clf = MLPClassifier(hidden_layer_sizes=tuple([neuron]*layer),
                                    max_iter=1000,
                                    activation="tanh",
                                    solver="adam",
                                    random_state=r,
                                    alpha=0.0001,
                                    early_stopping=True)
                clf.fit(data_train,target_train)
                score = clf.score(data_test,target_test)
                score_list.append(score)
                
            score_ave,score_std = np.average(score_list),np.std(score_list)
            print(neuron,layer," -> ",score_ave,score_std)
            
            if score_std/score_ave <= 0.1 and s_best < score_ave:
                n_best = neuron
                l_best = layer
                s_best = score_ave
        
        pen,neu,lay,act = 0.0001,n_best,l_best,"tanh"
        
 
        # bounds = [{"name":"penalty","type":"continuous","domain":(0.00001,0.0001)},
        #           {"name":"neurons","type":"discrete","domain":(3,4,5,6,7,8,9,10)},
        #           {"name":"layers","type":"discrete","domain":(1,2,3,4,5)},]
        # activations = ["identity","logistic","tanh","relu"]
        
        # def func_NN(x):
        #     penalty = float(x[:,0][0])
        #     neuron = int(x[:,1][0])
        #     layer = int(x[:,2][0])
        #     activation = activations[2]

        #     score_list = []
        #     for r in [0,1,2,3,4]:
        #         data_train,data_test,target_train,target_test = train_test_split(ans_input,ans_raw["label"],
        #                                                                          test_size=0.4,random_state=r,stratify=ans_raw["label"])
        #         clf = MLPClassifier(hidden_layer_sizes=tuple([neuron]*layer),
        #                             max_iter=1000,
        #                             activation=activation,
        #                             solver="adam",
        #                             random_state=r,
        #                             alpha=penalty,
        #                             early_stopping=True)
        #         clf.fit(data_train,target_train)
        #         score = clf.score(data_test,target_test)
        #         score_list.append(score)
            
        #     score_ave,score_std = np.average(score_list),np.std(score_list)
        #     print(neuron,layer,penalty,activation," -> ",score_ave,score_std)
        #     if score_std/score_ave > 0.1:
        #         return 0
        #     return -1 * score_ave
        
        # opt_NN = GPyOpt.methods.BayesianOptimization(f=func_NN,domain=bounds)
        # opt_NN.run_optimization(max_iter=100)
        # res = opt_NN.x_opt
        # pen,neu,lay,act = float(res[0]),int(res[1]),int(res[2]),activations[2]


        print("hyper parameters :",pen,neu,lay,act)

        f = open(args.savepath + "/" + args.foldername + "/" + "parameters.txt", "a")
        f.write('# of neurons : %s\n' %str(neu))
        f.write('# of hidden layers : %s\n' %str(lay))
        f.write('solver : adam\n')
        f.write('activation function : %s\n' %str(act))
        f.write('penalty : %s\n' %str(pen))
        f.close()

        #reclassification
        noise_class = [0]*len(noise_raw)

        linetype = ["solid","dashed","dashdot","dotted",(0,(10,10))]
        fig, ax = plt.subplots(figsize=(6,4))
        for r in [0,1,2,3,4]:
            clf = MLPClassifier(hidden_layer_sizes=tuple([neu]*lay),
                                max_iter=1000,
                                activation=act,
                                solver="adam",
                                random_state=r,
                                alpha=pen,
                                early_stopping=True)
            clf.fit(ans_input,ans_raw["label"])
            class_prob = clf.predict_proba(noise_input)
            max_prob_list = []
            for l in class_prob:
                max_prob_list.append(max(l))
            for i in range(len(max_prob_list)):
                if noise_class[i] != -1:
                    if max_prob_list[i] >= const.probability:
                        c = list(class_prob[i]).index(max_prob_list[i])+1
                        if noise_class[i] == 0 or noise_class[i] == c:
                            noise_class[i] = c
                        else:
                            noise_class[i] = -1
                    else:
                        noise_class[i] = -1
            ax.plot(clf.loss_curve_,linestyle=linetype[r],c="orangered",linewidth=1)
            ax.plot(clf.validation_scores_,linestyle=linetype[r],c="dodgerblue",linewidth=1)
        ax.set_xlabel('Number of iterations')
        ax.set_ylabel('Loss or accuracy')
        plt.savefig(self.args.savepath + "/" + self.args.foldername + "/clustering/learning_curve.png",bbox_inches="tight",pad_inches=0.05)
        plt.clf()
        plt.close()

        noise_raw["label"] = noise_class
        tremorB = noise_raw[noise_raw["label"]==tremor_clu_num].copy()
        
        if len(tremorB) == 0:
            print("no Quality B events")
        else:
            fig = md.hypo_colorbar(tremorB["lat"],tremorB["lon"],
                                tremorB["dep"],"Depth (km)",
                                lat_range=const.lat_range,lon_range=const.lon_range)
            fig.savefig(self.args.savepath + "/" + self.args.foldername + "/clustering/QualityB.png")

            tremorA["quality"] = "A"
            tremorB["quality"] = "B"
            
            tremor_catalog = pd.concat([tremorA,tremorB])
            tremor_catalog.sort_values(by ="ymdhms", ascending = True, inplace = True)
            tremor_catalog.reset_index(drop=True, inplace=True)
            columns = self.columns
            columns.append("quality")
            tremor_catalog = tremor_catalog[columns].copy()
            tremor_catalog.to_csv(self.args.savepath + "/" + self.args.foldername + "/tremor_catalog.csv",index=False)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--year', help='year or all',required=True)
    parser.add_argument('--savepath', help='save path',required=True)
    parser.add_argument('--foldername',help='folder name',required=True)
    args = parser.parse_args()

#    subprocess.call("rm %s/%s/clustering/*" %(args.savepath,args.foldername),shell=True)

    events(args)
    
    Cl = Clustering(args)
    Cl.clustering()
    Cl.compile_AB()
    Cl.reclassification()