# make_figure.py
# Simon Hulse
# simonhulse@protonmail.com
# Last Edited: Wed 24 Jan 2024 19:41:42 EST

import matplotlib as mpl
import matplotlib.pyplot as plt
import nmrespy as ne
import numpy as np


class JresAB(ne.ExpInfo):

    def __init__(
        self,
        pts1: int,
        pts2: int,
        sw1: float,
        sw2: float,
        offset: float,
        omega_a: float,
        omega_b: float,
        J_ab: float,
    ) -> None:
        super().__init__(dim=2, sw=(sw1, sw2), offset=(0., offset), default_pts=(pts1, pts2), nuclei=(None, "1H"))
        # N.B. omega_a is greater than omega_b
        assert omega_a != omega_b
        self.omega_a = max(omega_a, omega_b)
        self.omega_b = min(omega_a, omega_b)
        self.J_ab = J_ab
        self.half_J = 0.5 * self.J_ab
        self.theta = 0.5 * np.arctan(J_ab / (omega_a - omega_b))
        self.delta_omega = np.sqrt((self.omega_a - self.omega_b) ** 2 + self.J_ab ** 2)
        self.half_delta_omega = 0.5 * self.delta_omega

    @property
    def first_order_info(self):
        return [
            (0.5, (self.half_J, self.omega_a + self.half_J)),
            (0.5, (-self.half_J, self.omega_a - self.half_J)),
            (0.5, (self.half_J, self.omega_b + self.half_J)),
            (0.5, (-self.half_J, self.omega_b - self.half_J)),
        ]

    @property
    def artefact_info(self):
        return [
            (self.theta, (self.half_J + self.half_delta_omega, self.omega_a + self.half_J)),
            (-self.theta, (-self.half_J + self.half_delta_omega, self.omega_a - self.half_J)),
            (-self.theta, (self.half_J - self.half_delta_omega, self.omega_b + self.half_J)),
            (self.theta, (-self.half_J - self.half_delta_omega, self.omega_b - self.half_J)),
        ]

    def first_order_spectrum(self, lb: float = 10.0, shear: bool = False):
        return self._make_spectrum(self.first_order_info, lb, shear)

    def artefact_spectrum(self, lb: float = 10.0, shear: bool = False):
        return self._make_spectrum(self.artefact_info, lb, shear)

    def first_order_artefact_pair_spectra(self, lb: float = 10.0, shear: bool = False):
        for i, (finfo, ainfo) in enumerate(zip(self.first_order_info, self.artefact_info)):
            info = [finfo, ainfo]
            if i == 0:
                f1, f2, spec = self._make_spectrum(info, lb, shear)
                spectra = np.zeros((4, *spec.shape), dtype=spec.dtype)
            else:
                _, _, spec = self._make_spectrum(info, lb, shear)
            spectra[i] = spec
        return f1, f2, spectra

    def full_spectrum(self, lb: float = 10.0, shear: bool = False):
        f1, f2, first_order = self._make_spectrum(self.first_order_info, lb, shear)
        _, _, artefact = self._make_spectrum(self.artefact_info, lb, shear)
        spectrum = first_order + artefact
        return f1, f2, spectrum

    def _make_spectrum(self, info, lb: float, shear: bool):
        fid = np.zeros(self.default_pts, dtype="complex128")
        params = np.zeros((4, 6), dtype="float64")
        for i, (amp, (f1, f2)) in enumerate(info):
            params[i, 0] = amp
            params[i, 2] = f1
            params[i, 3] = f2
            if shear:
                params[i, 3] -= f1

        fid = self.make_fid(params)
        fid = ne.sig.exp_apodisation(fid, k=lb)
        fid = ne.sig.zf(fid)
        spectrum = np.abs(ne.sig.ft(ne.sig.sinebell_apodisation(fid)))
        f1, f2 = self.get_shifts(pts=spectrum.shape)
        return f1, f2, spectrum


pts1 = 128
pts2 = 128
sw1 = 200.0
sw2 = 220.0
offset = 0.0
J_ab = 20.0
diffs = [150, 100, 50]
fig, axs = plt.subplots(
    ncols=3,
    nrows=1,
    gridspec_kw={
        "left": 0.07,
        "right": 0.99,
        "bottom": 0.12,
        "wspace": 0.05,
    },
    figsize=(6., 2.5),
)
colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
base = 0.04
factor = 2.0
nlevels = 10
levels = [base * factor ** i for i in range(nlevels)]
textinfo = [
    dict(s="$f_\\text{A} + J_{\\text{AB}}$", y=-80.0, ha="left"),
    dict(s="$f_\\text{A} - J_{\\text{AB}}$", y=-60.0, ha="left"),
    dict(s="$f_\\text{B} + J_{\\text{AB}}$", y=60.0, ha="right"),
    dict(s="$f_\\text{B} - J_{\\text{AB}}$", y=80.0, ha="right"),
]
for i, (diff, ax) in enumerate(zip(diffs, axs)):
    omega_a, omega_b = offset + 0.5 * diff, offset - 0.5 * diff
    dataset = JresAB(pts1, pts2, sw1, sw2, offset, omega_a, omega_b, J_ab)
    f1, f2, spectrum = dataset.full_spectrum()
    ax.contour(f2, f1, spectrum, colors="k", levels=levels, linewidths=0.4)
    for j, (finfo, ainfo) in enumerate(zip(dataset.first_order_info, dataset.artefact_info)):
        f_f1, f_f2 = finfo[1]
        a_f1, a_f2 = ainfo[1]
        ax.axvline(f_f2, color=colors[j], zorder=-1)
        if i == 0:
            tinfo = textinfo[j]
            hpad = 2.5
            if tinfo["ha"] == "left":
                hpad = -hpad - 2.0
            ax.text(
                f_f2 + hpad,
                tinfo["y"],
                tinfo["s"],
                ha=tinfo["ha"],
                va="center",
                color=colors[j],
                fontsize=7,
                bbox=dict(
                    facecolor="w",
                    edgecolor="none",
                    pad=1.,
                )
            )
            if j == 1:
                xpos = 20.0
                ax.plot([a_f2, xpos], [a_f1, a_f1], color="k", ls=":")
                ax.plot([a_f2, xpos], [f_f1, f_f1], color="k", ls=":")
                arrow = mpl.patches.FancyArrowPatch([xpos, a_f1], [xpos, f_f1], arrowstyle="<->", shrinkA=0.0, shrinkB=0.0, color="k", mutation_scale=10.0, lw=0.8)
                ax.add_patch(arrow)
                ax.text(xpos - 2.0, 0.5 * (a_f1 + f_f1), "$\\xi_{\\text{AB}}$", va="center", fontsize=7)
    ax.set_xlim(offset + 0.5 * sw2, offset - 0.5 * sw2)
    ax.set_ylim(0.5 * sw1, -0.5 * sw1)
    if i != 0:
        ax.set_yticks([])
    else:
        ax.set_ylabel("$F^{(1)}$ (Hz)")
    if i == 1:
        ax.set_xlabel("$F^{(2)}$ (Hz)")
    ax.text(0.02, 0.94, f"\\textbf{{{chr(97 + i)}.}}", transform=ax.transAxes)
axs[0].text(0.5, 0.5, "s", fontsize=7, ha="center", va="center")
axs[0].text(0.5, 0.6, "f", fontsize=7, ha="center", va="center")

fig.savefig("figures/jres_ab/jres_ab.pdf")
