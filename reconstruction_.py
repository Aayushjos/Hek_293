import numpy as np
import matplotlib.pyplot as plt
from frft22d import frft22d

Ho = plt.imread('C:/Users/EEE/Pictures/Experiment/230524/Image-26.tif').astype(float)
Hr = plt.imread('C:/Users/EEE/Pictures/Experiment/230524/Image-16.tif').astype(float)


Ho = Ho - np.mean(Hr)


Ho_upd = np.pad(Ho, ((250, 250), (250, 250)), mode='constant')


M, N = Ho_upd.shape


zeta = np.arange(1, M+1)
eta = np.arange(1, N+1)
Zeta, Eta = np.meshgrid(zeta, eta)


pixel_pitch = 5.5e-6 / 20  # pixel pitch in meters
wavelength = 532e-9  # wavelength in meters


h1 = M * pixel_pitch
h2 = N * pixel_pitch


plt.figure()
plt.imshow(Ho, cmap='gray')
plt.colorbar()
plt.show()


angles = np.arange(0.79, 0.804, 0.001)


frft_o2_1 = np.zeros((M, N, len(angles)), dtype=np.complex128)


for i, angle in enumerate(angles):
    
    frft_o2_1[:, :, i] = frft22d(Ho_upd, [angle, angle])
    plt.imshow(np.log(np.abs(frft_o2_1[:, :, i])), cmap='gray')
    plt.title(f"The angle is {angle}")
    plt.show()
    plt.pause(1)


Y_max = np.max(np.log(np.abs(frft_o2_1)), axis=(0, 1))


plt.stem(angles, Y_max)
plt.xlabel('Angle')
plt.ylabel('Max Value')
plt.title('Maximum value vs Angle')
plt.show()
