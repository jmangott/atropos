import numpy as np
import matplotlib.pyplot as plt

plt.style.use("./scripts/notebooks/custom_style_boolean.mplstyle")

xpos = [1.025, 1.45, 1.531, 2.25, 2.5, 3, 3.25, 4.061, 4.5]
r2 = [
    0.06249999658962048,
    0.06249993817654081,
    0.03580575595163796,
    0.015476746806795005,
    0.0624999381765408,
    0.06249998967422543,
    0.06250189374148563,
    0.061886865211785394,
    0.06249998175514748,
]
coef2 = np.polyfit(xpos, r2, 1)
reg2 = np.poly1d(coef2)

r4 = [
    0.02767572112065949,
    0.026219852191997732,
    0.01951168230570425,
    0.014793954536987966,
    0.010586538852791344,
    0.017148191797838004,
    0.02204514515890582,
    0.04603083842319218,
    0.06249993817654081,
]
coef4 = np.polyfit(xpos, r4, 1)
reg4 = np.poly1d(coef4)

r8 = [
    0.00370577261843439,
    0.012777562053622349,
    0.0028898680592915824,
    0.006971756682516046,
    0.00772123710302406,
    0.011809388434920118,
    0.015153961350549784,
    0.02639736306813105,
    0.01538224661495153,
]
coef8 = np.polyfit(xpos, r8, 1)
reg8 = np.poly1d(coef8)

r16 = [
    0.0008335001126779373,
    0.0014916257642706918,
    0.00095222559883653,
    0.0007451089484385581,
    0.0006213575890701081,
    0.0014081178255637844,
    0.0006893686065158333,
    0.01633032165986715,
    0.008881641466096851,
]
coef16 = np.polyfit(xpos, r16, 1)
reg16 = np.poly1d(coef16)

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

fig, ax1 = plt.subplots(figsize=(5, 4.5))
ax1.plot(1.025, 0.06249999658962048, "o", fillstyle="none", color=colors[0])
ax1.plot(1.025, 0.02767572112065949, "o", fillstyle="none", color=colors[1])
ax1.plot(1.025, 0.00370577261843439, "o", fillstyle="none", color=colors[2])
ax1.plot(1.025, 0.0008335001126779373, "o", fillstyle="none", color=colors[3])

ax1.plot(1.45, 0.06249993817654081, "o", fillstyle="none", color=colors[0])
ax1.plot(1.45, 0.026219852191997732, "o", fillstyle="none", color=colors[1])
ax1.plot(1.45, 0.012777562053622349, "o", fillstyle="none", color=colors[2])
ax1.plot(1.45, 0.0014916257642706918, "o", fillstyle="none", color=colors[3])

ax1.plot(1.531, 0.03580575595163796, "o", fillstyle="none", color=colors[0])
ax1.plot(1.531, 0.01951168230570425, "o", fillstyle="none", color=colors[1])
ax1.plot(1.531, 0.0028898680592915824, "o", fillstyle="none", color=colors[2])
ax1.plot(1.531, 0.00095222559883653, "o", fillstyle="none", color=colors[3])

ax1.plot(2.25, 0.015476746806795005, "x", fillstyle="none", color=colors[0])
ax1.plot(2.25, 0.014793954536987966, "x", fillstyle="none", color=colors[1])
ax1.plot(2.25, 0.006971756682516046, "x", fillstyle="none", color=colors[2])
ax1.plot(2.25, 0.0007451089484385581, "x", fillstyle="none", color=colors[3])

ax1.plot(2.5, 0.0624999381765408, "x", fillstyle="none", color=colors[0])
ax1.plot(2.5, 0.010586538852791344, "x", fillstyle="none", color=colors[1])
ax1.plot(2.5, 0.00772123710302406, "x", fillstyle="none", color=colors[2])
ax1.plot(2.5, 0.0006213575890701081, "x", fillstyle="none", color=colors[3])

ax1.plot(3, 0.06249998967422543, "x", fillstyle="none", color=colors[0])
ax1.plot(3, 0.017148191797838004, "x", fillstyle="none", color=colors[1])
ax1.plot(3, 0.011809388434920118, "x", fillstyle="none", color=colors[2])
ax1.plot(3, 0.0014081178255637844, "x", fillstyle="none", color=colors[3])

ax1.plot(3.25, 0.06250189374148563, "o", fillstyle="none", color=colors[0])
ax1.plot(3.25, 0.02204514515890582, "o", fillstyle="none", color=colors[1])
ax1.plot(3.25, 0.015153961350549784, "o", fillstyle="none", color=colors[2])
ax1.plot(3.25, 0.0006893686065158333, "o", fillstyle="none", color=colors[3])

ax1.plot(4.061, 0.061886865211785394, "^", fillstyle="none", color=colors[0])
ax1.plot(4.061, 0.04603083842319218, "^", fillstyle="none", color=colors[1])
ax1.plot(4.061, 0.02639736306813105, "^", fillstyle="none", color=colors[2])
ax1.plot(4.061, 0.01633032165986715, "^", fillstyle="none", color=colors[3])

ax1.plot(4.5, 0.06249998175514748, "+", fillstyle="none", color=colors[0])
ax1.plot(4.5, 0.06249993817654081, "+", fillstyle="none", color=colors[1])
ax1.plot(4.5, 0.01538224661495153, "+", fillstyle="none", color=colors[2])
ax1.plot(4.5, 0.008881641466096851, "+", fillstyle="none", color=colors[3])


(l0,) = ax1.plot(xpos, reg2(xpos), color=colors[0], label="rank 2")
(l1,) = ax1.plot(xpos, reg4(xpos), color=colors[1], label="rank 4")
(l2,) = ax1.plot(xpos, reg8(xpos), color=colors[2], label="rank 8")
(l4,) = ax1.plot(xpos, reg16(xpos), color=colors[3], label="rank 16")

markerstyles = ["x", "o", "+", "^"]
markerlabels = ["3 cuts", "4 cuts", "5 cuts", "6 cuts"]

dummy_markers = []
for marker in markerstyles:
    dummy_markers.append(
        ax1.plot([], [], ls="none", marker=marker, fillstyle="none", color="gray")[0]
    )

linehandles, linelabels = ax1.get_legend_handles_labels()
legend1 = fig.legend(
    linehandles,
    linelabels,
    ncol=len(linelabels),
    loc="center",
    bbox_to_anchor=(0.5, 1.06),
)
legend2 = fig.legend(
    dummy_markers,
    markerlabels,
    ncol=len(markerlabels),
    loc="center",
    handlelength=0.3,
    bbox_to_anchor=(0.5, 1.01),
)
fig.add_artist(legend1)

ax1.set_xlabel("entropy $H$")

xpos_adjusted = xpos
xpos_adjusted[2] = xpos[2]+0.1

ax1.set_xticks(xpos_adjusted, xpos, rotation=90)

ax2 = ax1.twiny()
ax1.set_xlim([0.925, 4.6])
ax1.set_ylim([-0.005, 0.065])
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(
    xpos,
    ["M11", "M4", "PURPLE", "M1", "BLUE", "M3", "M10", "M32", "GREEN"],
    rotation=90,
)
ax2.set_xlabel("partitioning")

ax1.set_ylabel("$\| P - P_\mathrm{{ref}}\|_\infty$")
plt.tight_layout()
plt.savefig("plots/mTor_figure2.pdf", bbox_inches="tight")
