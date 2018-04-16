import sys
import numpy as numpy
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

cols = ['sensor_type', 'nis']
with open('NIS_data.csv') as f:
  nis_data = pd.read_csv(f, names=cols, lineterminator='\n')

nis_laser = []
nis_radar = []

for i in range(len(nis_data)):
  if nis_data['sensor_type'][i] == 'LASER':
    nis_laser.append(nis_data['nis'][i])
  else:
    nis_radar.append(nis_data['nis'][i])

confidence95Radar = 7.82
confidence95Lidar = 5.99

plt.figure(figsize=(8,4))
plt.plot(nis_laser)
plt.ylabel('NIS')
plt.xlabel('Time Step')
plt.title('Lidar NIS')
plt.axhline(y=confidence95Lidar, color='r', linestyle='-')
plt.ylim([0,15])
plt.tight_layout()
plt.savefig("NIS_Laser.png", bbox_inches='tight', dpi=300)

plt.figure(figsize=(8,4))
plt.plot(nis_radar)
plt.ylabel('NIS')
plt.xlabel('Time Step')
plt.title('Radar NIS')
plt.axhline(y=confidence95Radar, color='r', linestyle='-')
plt.ylim([0,15])
plt.tight_layout()
plt.savefig("NIS_Radar.png", bbox_inches='tight', dpi=300)