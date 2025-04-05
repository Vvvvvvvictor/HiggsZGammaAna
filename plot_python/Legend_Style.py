from matplotlib import font_manager
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(6, 6))

# Dummy plot
ax.plot([0, 1], [0, 1], label="Serif")
ax.plot([0, 1], [1, 0], label="Sans-serif")
ax.plot([0, 1], [0.5, 0.5], label="Monospace")
ax.plot([0.5, 0.5], [0, 1], label="Cursive")
ax.plot([0.2, 0.8], [0.2, 0.8], label="Fantasy")

# Define different font styles
serif_font = font_manager.FontProperties(family="serif", size=14, style="italic")
sans_font = font_manager.FontProperties(family="sans-serif", size=14, weight="bold")
mono_font = font_manager.FontProperties(family="monospace", size=14)
cursive_font = font_manager.FontProperties(family="cursive", size=14, style="italic")
fantasy_font = font_manager.FontProperties(family="fantasy", size=14, weight="bold")

# Apply multiple legends at different positions
leg1 = ax.legend(loc="upper left", prop=serif_font, title="Serif Font")
leg2 = ax.legend(loc="upper right", prop=sans_font, title="Sans-serif Font")
leg3 = ax.legend(loc="lower left", prop=mono_font, title="Monospace Font")
leg4 = ax.legend(loc="lower right", prop=cursive_font, title="Cursive Font")
leg5 = ax.legend(loc="center", prop=fantasy_font, title="Fantasy Font")

# Add all legends back to the figure
ax.add_artist(leg1)
ax.add_artist(leg2)
ax.add_artist(leg3)
ax.add_artist(leg4)

# Save and show the figure
plt.savefig("pic/all_legend_styles.pdf", format="pdf", dpi=300, bbox_inches="tight")
