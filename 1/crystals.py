import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import os
import glob

def create_circle(cx, cy, r=0.5, resolution=20):
    theta = np.linspace(0, 2 * np.pi, resolution)
    x = r * np.cos(theta) + cx
    y = r * np.sin(theta) + cy
    return x, y

# Define input and output directories
input_dir = "output/crystals"
output_dir = "output/plotted_crystals"
os.makedirs(output_dir, exist_ok=True)

# Get all .txt files in the input directory
file_list = glob.glob(os.path.join(input_dir, "*.txt"))

for file_path in file_list:
    try:
        # Extract filename and density
        file_name = os.path.basename(file_path)
        density = float(file_name.replace('.txt', ''))
        
        # Define box size (a) based on density
        a = (40 / density) ** 0.5
        d = 1  # Circle diameter

        # Read the coordinates from the file
        data = []
        with open(file_path, "r") as file:
            for line in file:
                if line.strip() and not line.lower().startswith("packing"):
                    values = line.split()
                    if len(values) == 2:  # Only process lines with two coordinates (2D)
                        data.append([float(values[0]), float(values[1])])

        # Convert list to numpy array
        coords = np.array(data)

        # Check for overlap
        overlap_found = False
        for i in range(len(coords)):
            for j in range(i + 1, len(coords)):
                if np.linalg.norm(coords[i] - coords[j]) < d:
                    print(f"Overlap detected in {file_name} between {coords[i]} and {coords[j]}")
                    overlap_found = True
        
        if not overlap_found:
            print(f"No overlap detected in {file_name}.")

        # Create figure
        fig = go.Figure()

        # Add circles
        for cx, cy in coords:
            if abs(cx) <= a/2 and abs(cy) <= a/2:
                x, y = create_circle(cx, cy, r=d/2)
                fig.add_trace(go.Scatter(
                    x=x, y=y,
                    fill="toself",
                    mode="lines",
                    line=dict(color="blue"),
                    opacity=0.8,
                    showlegend=False
                ))

        fig.update_layout(
            xaxis=dict(range=[-a/2, a/2], fixedrange=True, showgrid=True, zeroline=False),
            yaxis=dict(range=[-a/2, a/2], fixedrange=True, showgrid=True, zeroline=False),
            title=f"Packing: 2D Circle Packing (Density {density})",
            width=600,
            height=600,
            plot_bgcolor="white",
            paper_bgcolor="white",
        )

        # Save plot as an image in the output directory
        output_file = os.path.join(output_dir, file_name.replace('.txt', '.png'))
        pio.write_image(fig, output_file)
        print(f"Saved plot: {output_file}")
    
    except Exception as e:
        print(f"Error processing {file_name}: {e}")
