import matplotlib.pyplot as plb
import pandas as pd
from osm import OSM as OSM
from fps_style import *

cat = pd.read_csv('eq_until_2M_text.csv')
fig1 = plb.figure()
ax1 = fig1.add_subplot(1, 2, 1)
ax2 = fig1.add_subplot(1, 2, 2)
dP, dT, dB = convert2tri(cat.Strike.values, cat.Dip.values, cat.Rake.values)
mClass = TernaryFMS(dP, dT, dB, ax1)

dxy = 0.2
ax2.set_xlim([min(cat.Longitude.values)-dxy, max(cat.Longitude.values)+dxy])
ax2.set_ylim([min(cat.Latitude.values)-dxy, max(cat.Latitude.values)+dxy])

cm = ['grey', 'red', 'green', 'blue']
labels = ['Oblique', 'Strike-slip', 'Normal', 'Reverse']
for ii in range(len(labels)):
    ax2.scatter(0, 0, c=cm[ii], label=labels[ii])
for ii in range(len(dP)):
    b2 = beach([cat.Strike.values[ii], cat.Dip.values[ii], cat.Rake.values[ii]], xy=(cat.Longitude.values[ii], cat.Latitude.values[ii]), width=0.005 * cat.Magnitude.values[ii]**1.5, linewidth=0.2, facecolor=cm[int(mClass[ii])])
    ax2.add_collection(b2)
ax2.legend()
ax2.osm = OSM(ax2)
ax2.osm.draw()
plb.show()