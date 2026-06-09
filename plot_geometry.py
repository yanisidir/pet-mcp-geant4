import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Dimensions en cm
mcp_radius = 2.5      # 25 mm
mcp_length = 0.3      # 3 mm
mcp_gap = 0.005       # 50 um = 0.005 cm
first_mcp_z = 20.0
n_mcp = 2

fig, ax = plt.subplots(figsize=(10, 4))

# Source
ax.scatter(0, 0, s=80, label="Source", marker="*")

# MCP stacks
for side, color, label in [(1, "tab:blue", "+z stack"), (-1, "tab:green", "-z stack")]:
    for i in range(n_mcp):
        z = side * (first_mcp_z + i * (mcp_length + mcp_gap))

        rect = Rectangle(
            (z - side * mcp_length / 2, -mcp_radius),
            side * mcp_length,
            2 * mcp_radius,
            alpha=0.5,
            label=label if i == 0 else None
        )

        ax.add_patch(rect)
        ax.text(z, mcp_radius + 0.3, f"MCP {i}", ha="center", fontsize=10)

# Photons back-to-back
ax.arrow(0, 0, 18, 0, head_width=0.25, head_length=0.6, length_includes_head=True)
ax.arrow(0, 0, -18, 0, head_width=0.25, head_length=0.6, length_includes_head=True)
ax.text(8, 0.3, "511 keV gamma", ha="center")
ax.text(-8, 0.3, "511 keV gamma", ha="center")

ax.set_xlabel("z position [cm]")
ax.set_ylabel("Transverse size [cm]")
ax.set_title("PET-like Geant4 geometry: two opposite MCP stacks")
ax.set_aspect("equal", adjustable="box")
ax.grid(True)
ax.legend()
ax.set_xlim(-22, 22)
ax.set_ylim(-4, 4)

plt.tight_layout()
plt.savefig("mcp_pet_geometry_side.png", dpi=300)
plt.savefig("mcp_pet_geometry_side.pdf")
plt.show()