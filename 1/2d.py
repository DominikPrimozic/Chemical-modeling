import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
# Define box size (a) and sphere diameter
a = (40/0.5)**(1/3)  # Adjust as needed
d = 1  # Sphere diameter

# Read the coordinates from the file
data = []
with open("output/placements/3.40.0.500000.txt", "r") as file:
    for line in file:
        if line.strip() and not line.lower().startswith("packing"):
            values = line.split()
            if len(values) == 3:  # Only process lines with three coordinates
                data.append([float(values[0]), float(values[1]), float(values[2])])

# Convert list to numpy array
coords = np.array(data)

# Create sphere mesh
def create_sphere(cx, cy, cz, r=0.5, resolution=20):
    u = np.linspace(0, np.pi, resolution)
    v = np.linspace(0, 2 * np.pi, resolution)
    x = r * np.outer(np.sin(u), np.cos(v)) + cx
    y = r * np.outer(np.sin(u), np.sin(v)) + cy
    z = r * np.outer(np.cos(u), np.ones_like(v)) + cz
    return x, y, z

# Create figure
fig = go.Figure()

# Add spheres
for cx, cy, cz in coords:
    x, y, z = create_sphere(cx, cy, cz, r=d/2)
    fig.add_trace(go.Surface(
        x=x, y=y, z=z,
        colorscale=[[0, 'blue'], [1, 'blue']],
        showscale=False,
        opacity=0.8
    ))

# Set layout
fig.update_layout(
    scene=dict(
        xaxis=dict(range=[-a, a]),
        yaxis=dict(range=[-a, a]),
        zaxis=dict(range=[-a, a]),
        aspectmode='cube'
    )
)

# Open in browser
pio.show(fig, renderer='browser')

# Function to calculate distance between two points
def distance(p1, p2):
    return np.linalg.norm(p1 - p2)

# Check for overlap
overlap_found = False
for i in range(len(coords)):
    for j in range(i + 1, len(coords)):  # Compare each pair only once
        dist = distance(coords[i], coords[j])
        if dist < d:  # If distance is less than diameter, there is overlap
            print(f"Overlap detected between spheres at {coords[i]} and {coords[j]}")
            overlap_found = True

if not overlap_found:
    print("No overlap detected.")