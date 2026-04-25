import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ---------------- USER SETTINGS ----------------
file_name = "TimeStep.out"
N = 128
interval_ms = 150

# ---------------- READ TIMESTEPS ----------------
def read_frames(file_name, N):
    frames = []
    current = []

    with open(file_name, "r") as f:
        for line in f:
            line = line.strip()

            if line == "":
                if len(current) == N:
                    frames.append(np.array(current, dtype=float))
                current = []
                continue

            parts = line.split()

            try:
                row = [float(v) for v in parts]
            except ValueError:
                continue

            if len(row) == N:
                current.append(row)

        if len(current) == N:
            frames.append(np.array(current, dtype=float))

    return frames

frames = read_frames(file_name, N)

if len(frames) == 0:
    raise RuntimeError("No 64x64 frames found. Check file format or N.")

print(f"Loaded {len(frames)} frames.")

# ---------------- ANIMATION ----------------
fig, ax = plt.subplots()

vmin = min(frame.min() for frame in frames)
vmax = max(frame.max() for frame in frames)

im = ax.imshow(
    frames[0],
    origin="lower",
    extent=[0, 1, 0, 1],
    vmin=vmin,
    vmax=vmax,
    animated=True
)

cbar = fig.colorbar(im, ax=ax)
cbar.set_label("Field value")

ax.set_xlabel("x")
ax.set_ylabel("y")
title = ax.set_title("Timestep 0")

def update(frame_id):
    im.set_array(frames[frame_id])
    title.set_text(f"Timestep {frame_id}")
    return [im, title]

ani = FuncAnimation(
    fig,
    update,
    frames=len(frames),
    interval=interval_ms,
    blit=False,
    repeat=True
)

plt.show()