from pymol.cgo import CYLINDER
from pymol import cmd

def spectrumbar_rmsd(name="rmsd_bar", start="0,0,0", length=15.0, radius=0.3,
                     steps=5, colors=("red", "blue"), rmsd_min=10.0, rmsd_max=95.0):
    """
    Draws a labeled RMSD color gradient bar in PyMOL using CGO.

    Parameters:
        name:     Name of the CGO object
        start:    "x,y,z" string (use quotes in PyMOL)
        length:   Length of the gradient bar
        radius:   Radius of the cylinder
        steps:    Number of label intervals
        colors:   Tuple of PyMOL color names (low to high RMSD)
        rmsd_min: Minimum RMSD value (left)
        rmsd_max: Maximum RMSD value (right)
    """
    # Safely parse numeric parameters
    length = float(length)
    radius = float(radius)
    steps = int(steps)
    rmsd_min = float(rmsd_min)
    rmsd_max = float(rmsd_max)

    # Parse position string
    if isinstance(start, str):
        start = start.strip("\"'() ")
        x0, y0, z0 = map(float, start.split(","))
    else:
        x0, y0, z0 = start

    segments = 50
    dx = length / segments

    rgb_start = cmd.get_color_tuple(colors[0])
    rgb_end = cmd.get_color_tuple(colors[1])

    bar = []

    for i in range(segments):
        frac = i / segments
        r = rgb_start[0] * (1 - frac) + rgb_end[0] * frac
        g = rgb_start[1] * (1 - frac) + rgb_end[1] * frac
        b = rgb_start[2] * (1 - frac) + rgb_end[2] * frac

        x1 = x0 + i * dx
        x2 = x0 + (i + 1) * dx

        bar.extend([
            CYLINDER,
            x1, y0, z0,
            x2, y0, z0,
            radius,
            r, g, b,
            r, g, b
        ])

    cmd.delete(name)
    cmd.load_cgo(bar, name)

    # Add labels below the bar
    for j in range(steps + 1):
        frac = j / steps
        xpos = x0 + frac * length
        value = rmsd_min + frac * (rmsd_max - rmsd_min)
        cmd.pseudoatom(object=f"{name}_label_{j}", pos=[xpos, y0 - 1.0, z0], label=f"{value:.1f} Å")

    # Style labels
    cmd.set("label_size", 14)
    cmd.set("label_font_id", 7)

# Register with PyMOL
cmd.extend("spectrumbar_rmsd", spectrumbar_rmsd)
