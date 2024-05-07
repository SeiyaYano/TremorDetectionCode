import numpy as np
import modules as md
import math
import pandas as pd
from sklearn.linear_model import LinearRegression

def box_counting(df):
    
    print("calculating fractal dimension ...")
    df.reset_index(drop=True,inplace=True)
    
    origin = [np.average(df["lat"]), np.average(df["lon"])]
    xyz, xyz_list = [], []
    for i in range(len(df)):
        res = md.latlon2xy(origin,[df["lat"][i],df["lon"][i]])
        c = [res["x"],res["y"],abs(df["dep"][i])*-1]
        xyz.append(c)
        c_floor = [math.floor(res["x"]),math.floor(res["y"]),math.floor(abs(df["dep"][i])*-1)]
        if c_floor not in xyz_list:
            xyz_list.append(c_floor)
    
    df_xyz = pd.DataFrame(xyz)

    xmax, xmin = math.ceil(max(df_xyz[0]))+1, math.floor(min(df_xyz[0]))
    ymax, ymin = math.ceil(max(df_xyz[1]))+1, math.floor(min(df_xyz[1]))
    zmax, zmin = math.ceil(max(df_xyz[2]))+1, math.floor(min(df_xyz[2]))
    cube_len = max(xmax-xmin,ymax-ymin,zmax-zmin)

    array3d = np.zeros((cube_len, cube_len, cube_len))
    for c_floor in xyz_list:
        array3d[c_floor[0],c_floor[1],c_floor[2]] = 1

    def count(arr,x,y,z,b):
        i = j = k = 0
        b_tmp = b
        b_tmp2 = b
        ct = 0

        j_tmp = b
        k_tmp = b

        while k < z:
            while j < y:
                while i < x:
                    if (np.any(arr[i:b,j:j_tmp,k:k_tmp] == 1)):
                        ct += 1
                    i += b_tmp
                    b += b_tmp
                j += b_tmp2
                j_tmp += b_tmp2
                b = b_tmp2
                i = 0
            k += b_tmp2
            k_tmp += b_tmp2
            b = b_tmp2
            j = b_tmp2
            i =j = 0
        return ct

    graph_x, graph_y = [], []
    ct_1 = np.count_nonzero(array3d==1)
    grid_size = cube_len

    while grid_size >= 4:
        n = count(array3d,cube_len,cube_len,cube_len,grid_size)
        graph_x.append(math.log(grid_size))
        graph_y.append(math.log(n))
        grid_size = int(grid_size / 2)

    graph_x = np.array(graph_x).reshape((len(graph_x), 1))
    graph_y = np.array(graph_y).reshape((len(graph_y), 1))

    model_lr = LinearRegression()
    model_lr.fit(graph_x, graph_y)
    fractal = model_lr.coef_[0][0] * -1
    
    return fractal

