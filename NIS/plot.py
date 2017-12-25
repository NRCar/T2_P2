import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Normalized innovation squared (NIS) data values for radar and laser measurements but skip first NIS 
dataRadar = np.loadtxt( "radar_NIS.txt", skiprows=1)
dataLaser = np.loadtxt( "laser_NIS.txt",skiprows=1)

nisRadar = np.transpose(dataRadar)
nisLaser = np.transpose(dataLaser)
nis = [nisRadar, nisLaser]

# Laser measurements have 3 degrees of freedom.
# confidence95Laser is the threshold at which the p-value of the 
# NIS distribution is 0.05 i.e when we plot the NIS values for the
# set of radar measurements, roughly 5% of them should be
# above the 7.82 threshold, if the noise values are consistent.
# 
# This serves as a check on our choice of process noise values.
confidence95Radar = 7.82

# Laser measurements have 2 degrees of freedom,
# so the threshold is different.
confidence95Laser = 5.99

confidences = np.array([confidence95Radar, confidence95Laser])

fig = plt.figure(figsize=(16,8))
axes = []

id = 1
for data, confidence in zip(nis, confidences):
    subplot = 120 + id
    ax = fig.add_subplot(subplot)
    axes.append( ax )
    ax.tick_params(which='both',direction='in')
    ax.set_xlabel( "Measurement index")
    ax.plot( np.arange(0,len(data)), data, label="NIS",color='r')
    ax.axhline( y=confidence, color='b', linestyle='-'
        , label="confidence threshold")
    ax.legend()
    id += 1

axes[0].set_title("RADAR")
axes[1].set_title("LASER")

plt.tight_layout()
plt.savefig( "NIS.png", bbox_inches = 'tight', dpi = 300 )

plt.show()
