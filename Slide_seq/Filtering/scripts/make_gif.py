import os 
import imageio

png_folder = "/scratch/mfafouti/Mommybrain-PPD/Slide_seq/RCTD/Filtering/OUT/gif"
gif_path_base = "/scratch/mfafouti/Mommybrain-PPD/Slide_seq/RCTD/Filtering/OUT/looping_gifs"

images = sorted([f for f in os.listdir(png_folder) if f.endswith(".png")])
images = [os.path.join(png_folder, f) for f in images]  # full paths

gif_path = os.path.join(gif_path_base, "umap_animation.gif")
with imageio.get_writer(gif_path, mode='I', duration=300, loop=0) as writer:
    for img in images:
        frame = imageio.imread(img)
        writer.append_data(frame)

print(f"✅ GIF saved at {gif_path}")