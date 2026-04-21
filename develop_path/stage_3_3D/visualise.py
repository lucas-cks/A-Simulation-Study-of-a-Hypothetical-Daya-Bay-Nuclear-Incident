import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

NX, NY, NZ = 200, 200, 50
LX, LY, LZ = 100000.0, 100000.0, 5000.0   # m
DX, DY, DZ = LX/(NX-1), LY/(NY-1), LZ/(NZ-1)

x_km = np.linspace(0, LX/1000, NX)
y_km = np.linspace(0, LY/1000, NY)
z_km = np.linspace(0, LZ/1000, NZ)

filename = "output_3d.bin"
if not os.path.exists(filename):
    print(f"ERROR: {filename} not found.")
    exit(1)

frame_size = NX * NY * NZ * 8
file_size = os.path.getsize(filename)
if file_size % frame_size != 0:
    print(f"ERROR: File size {file_size} not multiple of frame size {frame_size}.")
    exit(1)

n_frames = file_size // frame_size
print(f"Number of frames: {n_frames}")

data = np.fromfile(filename, dtype=np.float64)
data = data.reshape(n_frames, NZ, NY, NX)
C = data[-1, :, :, :]   

# figure 1: ground concentration (z=0)
plt.figure(figsize=(10,8))
ground = C[0, :, :]
plt.imshow(ground, origin='lower', extent=[0, LX/1000, 0, LY/1000],
           cmap='jet', aspect='auto')
plt.colorbar(label='Concentration')
plt.xlabel('x (km)')
plt.ylabel('y (km)')
plt.title('Ground concentration (z=0)')
plt.savefig('ground_3d.png', dpi=150)
plt.show()

# figure 2: vertical slice at y=50 km
j_center = NY // 2
profile = C[:, j_center, :]
plt.figure(figsize=(12,5))
plt.imshow(profile, origin='lower', extent=[0, LX/1000, 0, LZ/1000],
           cmap='jet', aspect='auto')
plt.colorbar(label='Concentration')
plt.xlabel('x (km)')
plt.ylabel('z (km)')
plt.title(f'Vertical slice at y = {y_km[j_center]:.1f} km')
plt.savefig('vertical_3d.png', dpi=150)
plt.show()

# figure 3: 3D scatter of points above threshold
threshold = 0.05 * C.max()   
indices = np.where(C > threshold)
if len(indices[0]) > 0:
    step = 3   
    z_idx = indices[0][::step]
    y_idx = indices[1][::step]
    x_idx = indices[2][::step]
    vals = C[z_idx, y_idx, x_idx]
    xs = x_idx * DX / 1000
    ys = y_idx * DY / 1000
    zs = z_idx * DZ / 1000
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111, projection='3d')
    sc = ax.scatter(xs, ys, zs, c=vals, cmap='jet', s=3, alpha=0.5)
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_zlabel('z (km)')
    ax.set_title(f'3D plume (conc > {threshold:.2e})')
    plt.colorbar(sc, ax=ax, shrink=0.5)
    plt.savefig('scatter_3d.png', dpi=150)
    plt.show()
else:
    print("No points above threshold.")