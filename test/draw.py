"""
Allows to interactively draw a new polygon based
on an existing polygon.
"""
import sys
import matplotlib.pyplot as plt

in_file = sys.argv[-2]
out_file = sys.argv[-1]

fig = plt.figure()
fig = plt.figure()
ax = fig.add_subplot(111)

x = []
y = []
with open(in_file, 'r') as f:
    for line in f:
        x.append(float(line.split()[0]))    
        y.append(float(line.split()[1]))    

ax.plot(x, y)

def onclick(event):
    with open(out_file, 'a') as f:
        f.write("{0} {1}\n".format(event.xdata, event.ydata))
    print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          (event.button, event.x, event.y, event.xdata, event.ydata))

cid = fig.canvas.mpl_connect('button_press_event', onclick)

plt.show()
