import imageio
import os

def get_filenames(folder):
    names = []
    for currentdir, dirs, files in os.walk(folder):
        for filename in files:
                names.append(os.path.join(currentdir, filename))
    return names

folder = 'C:/TFZ/map/pictures'

pictures = get_filenames(folder)

with imageio.get_writer("./map/anime/anime_all.gif", fps=6) as writer:
    for picture in pictures:
        image = imageio.imread(picture)
        writer.append_data(image)
print("Gif saved")