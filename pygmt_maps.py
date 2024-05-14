import os
import pygmt
import pandas as pd
import numpy as np


def get_filenames(folder):
    names = []
    for currentdir, dirs, files in os.walk(folder):
        for filename in files:
                names.append(os.path.join(currentdir, filename))
    return names

# folder where all map-files are
mapfolder = './map/'

# folder where to save pictures
savefolder = './pictures/'

maps = get_filenames(mapfolder)

d = {}
        
for _map in maps:
    # gets date from the filename. Unnecessary to the program
    date = _map.split("/")[-1].split(".")[2] + '-' + _map.split("/")[-1].split(".")[6]
    fname = savefolder + date + ".png"
    print(fname)
    
    d["lon"] = []
    d["lat"] = []
    d["g"] = []
    with open(_map) as f:
        for line in f:
            lat, lon, g = map(float, line.split())
            d["lon"].append(180 - np.degrees(lon))
            d["lat"].append(90 - np.degrees(lat))
            d["g"].append(g)
            
            
    df = pd.DataFrame(data=d)
    
    region = [
            -180, 180,
            df.lat.min() - 1,
            df.lat.max() + 1
        ]
    

    pygmt.makecpt(cmap='plasma', series=[(df.g).min(), (df.g).max()])
    fig = pygmt.Figure()
    fig.basemap(region=region, projection="J-180/12c", frame=["afg", "+t"+date])
    fig.coast(land="black", water="skyblue")
    fig.plot(x=df.lon, y=df.lat, style='c0.035', fill=df.g, cmap=True)
    fig.colorbar(frame=['a', 'x+lg'])
    
    fig.savefig(fname)