# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Wed 19 Apr 2023 19:44:20 BST

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D, proj3d
from matplotlib.legend_handler import HandlerPatch

colours = plt.rcParams['axes.prop_cycle'].by_key()['color']

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))

        return np.min(zs)


fig = plt.figure(figsize = (6, 2.5))

left = 0
right = 0.99
xratio = 0.8
bottom = 0.11
top = 0.99
xgap = 0
ygap = 0.05
lwidth = (right - left - xgap) / (1 + 1 / xratio)
rwidth = lwidth / xratio
rheight = (top - bottom - ygap) / 2

hpad = 0.25
vpad = 0.22
axes = [
	fig.add_axes([0 - 1.05 * hpad, 0 - vpad, lwidth + 2 * hpad, 1  + 2 * vpad], projection='3d', computed_zorder=False),
	fig.add_axes([left+lwidth+xgap, bottom+rheight+ygap, rwidth, rheight]),
	fig.add_axes([left+lwidth+xgap, bottom, rwidth, rheight]),
]


time = 20
points = 2048
tp = np.linspace(0, time, points)

axes[0].plot(
    2 * np.sin(tp),
    2 * np.cos(tp),
    np.zeros(points),
    color="k",
    zorder=0,
)

coords = [
	([0,0], [0,0], [0, 1]),
	([0,0], [0,2.5], [0,0]),
	([0,2.5], [0,0], [0,0]),
	([0,0], [0,-2.5], [0,0]),
	([0,-2.5], [0,0], [0,0]),
]

texts = ['$z$', '$y$', '$x$', '$-y$', '$-x$']

for coord, text in zip(coords, texts):
	arrow = Arrow3D(
		*coord,
		mutation_scale=5,
		lw=1.5,
		arrowstyle="-|>",
		color="k",
		linestyle='-',
		zorder=0 if text != "$z$" else 600,
        shrinkA=0,
	)
	axes[0].add_artist(arrow)
	x, y, z = [c[1] for c in coord]
	if text == '$-x$':
		x -= 0.15
		y -= 0.4
	if text == '$-y$':
		x -= 0.15
		y -= 0.3

	axes[0].text(x, y, z, text)

a, f, d = 2., 2.7, 0.2
x = a * np.cos(f * tp) * np.exp(-tp * d)
y = - a * np.sin(f * tp) * np.exp(-tp * d)
z = 1 - np.exp(-tp * 0.1)

rate = 64
axes[0].plot(x, y, z, lw=1., zorder=700)
axes[0].scatter(
    x[::rate],
    y[::rate],
    z[::rate],
    zorder=500,
    depthshade=False,
    s=8,
)
axes[1].plot(tp, x)
axes[1].scatter(tp[::rate], x[::rate], s=8, edgecolor="k", lw=0.4, zorder=10)
axes[2].plot(tp, y)
axes[2].scatter(tp[::rate], y[::rate], s=8, edgecolor="k", lw=0.4, zorder=10)

# # Hide grid lines
axes[0].axis('off')

# Hide axes ticks
axes[0].set_xticks([])
axes[0].set_yticks([])
axes[0].set_zticks([])

axes[0].set_xlim(-2, 2)
axes[0].set_ylim(-2, 2)


axes[0].view_init(20, -45)
for ax in axes[1:]:
    ax.set_xticks(tp[::4 * rate])
    ax.set_ylim(-2.2, 2.2)
    ax.set_yticks([])
    ax.set_xlim(left=-1.2)
    xlim = ax.get_xlim()
    ax.plot(xlim, [0, 0], ls=":", color="k", zorder=0, lw=0.7)
    ax.set_xlim(xlim)

axes[2].set_xticklabels(
    [
        f"${i}\\Updelta_t$" if i != 0 else "$0$"
        for i in range(0, 4 * len(axes[2].get_xticks()), 4)
    ]
)
axes[1].set_xticklabels([])

axes[1].set_ylabel('$\\symbf{p}(x)$')
axes[2].set_ylabel('$\\symbf{p}(y)$')
axes[2].set_xlabel('$t$', labelpad=2)

for i, (x, y) in enumerate(zip((0.06, 0.445, 0.445), (0.94, 0.94, 0.48))):
    fig.text(x, y, f"\\textbf{{{chr(97 + i)}.}}")


fig.savefig("figures/quadrature_detection/quadrature_detection_new.pdf")
