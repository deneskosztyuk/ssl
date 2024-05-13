import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric, CartesianRepresentation

def plot_solar_system(t):
    # Set the reference frame to ICRS (International Celestial Reference System)
    with solar_system_ephemeris.set('builtin'):
        sun_pos = get_body_barycentric('sun', t)
        mercury_pos = get_body_barycentric('mercury', t)
        venus_pos = get_body_barycentric('venus', t)
        earth_pos = get_body_barycentric('earth', t)
        mars_pos = get_body_barycentric('mars', t)
        jupiter_pos = get_body_barycentric('jupiter', t)
        saturn_pos = get_body_barycentric('saturn', t)
        uranus_pos = get_body_barycentric('uranus', t)
        neptune_pos = get_body_barycentric('neptune', t)

    # Calculate relative positions
    mercury_rel_pos = (mercury_pos - sun_pos).represent_as(CartesianRepresentation).xyz.to(u.au)
    venus_rel_pos = (venus_pos - sun_pos).represent_as(CartesianRepresentation).xyz.to(u.au)
    earth_rel_pos = (earth_pos - sun_pos).represent_as(CartesianRepresentation).xyz.to(u.au)
    mars_rel_pos = (mars_pos - sun_pos).represent_as(CartesianRepresentation).xyz.to(u.au)
    jupiter_rel_pos = (jupiter_pos - sun_pos).represent_as(CartesianRepresentation).xyz.to(u.au)
    saturn_rel_pos = (saturn_pos - sun_pos).represent_as(CartesianRepresentation).xyz.to(u.au)
    uranus_rel_pos = (uranus_pos - sun_pos).represent_as(CartesianRepresentation).xyz.to(u.au)
    neptune_rel_pos = (neptune_pos - sun_pos).represent_as(CartesianRepresentation).xyz.to(u.au)

    # Clear the previous plot
    ax.clear()

    # Set up the plot
    ax.set_aspect('equal')
    ax.set_xlim(-10, 10)  # Adjust the limits to zoom in
    ax.set_ylim(-10, 10)  # Adjust the limits to zoom in
    ax.set_title('Solar System')
    ax.set_xlabel('X (au)')
    ax.set_ylabel('Y (au)')

    # Plot the Sun
    sun_radius = 0.2
    sun_color = 'yellow'
    sun = plt.Circle((0, 0), sun_radius, color=sun_color)
    ax.add_artist(sun)
    ax.annotate('Sun', xy=(0, 0), xytext=(10, 0), textcoords='offset points')

    # Plot the planets and their trajectories
    planets = [
        {'name': 'Mercury', 'pos': mercury_rel_pos, 'radius': 0.05, 'color': 'gray', 'trajectory_color': 'lightgray'},
        {'name': 'Venus', 'pos': venus_rel_pos, 'radius': 0.1, 'color': 'orange', 'trajectory_color': 'peachpuff'},
        {'name': 'Earth', 'pos': earth_rel_pos, 'radius': 0.1, 'color': 'blue', 'trajectory_color': 'skyblue'},
        {'name': 'Mars', 'pos': mars_rel_pos, 'radius': 0.08, 'color': 'red', 'trajectory_color': 'salmon'},
        {'name': 'Jupiter', 'pos': jupiter_rel_pos, 'radius': 0.2, 'color': 'brown', 'trajectory_color': 'sandybrown'},
        {'name': 'Saturn', 'pos': saturn_rel_pos, 'radius': 0.15, 'color': 'olive', 'trajectory_color': 'khaki'},
        {'name': 'Uranus', 'pos': uranus_rel_pos, 'radius': 0.12, 'color': 'lightblue', 'trajectory_color': 'powderblue'},
        {'name': 'Neptune', 'pos': neptune_rel_pos, 'radius': 0.12, 'color': 'darkblue', 'trajectory_color': 'lightskyblue'}
    ]

    for planet in planets:
        x, y = planet['pos'][0].value, planet['pos'][1].value
        planet_circle = plt.Circle((x, y), planet['radius'], color=planet['color'])
        ax.add_artist(planet_circle)
        ax.annotate(planet['name'], xy=(x, y), xytext=(10, 0), textcoords='offset points')

        # Plot the trajectory line
        trajectory_radius = np.sqrt(x**2 + y**2)
        angle = np.linspace(0, 2*np.pi, 100)
        trajectory_x = trajectory_radius * np.cos(angle)
        trajectory_y = trajectory_radius * np.sin(angle)
        ax.plot(trajectory_x, trajectory_y, linestyle='--', color=planet['trajectory_color'], linewidth=0.5)

    # Update the plot
    fig.canvas.draw_idle()

def animate(event):
    global t
    for _ in range(365):
        t += TimeDelta(1, format='jd')
        plot_solar_system(t)
        plt.pause(0.01)

# Set the initial time
t = Time("2023-05-14 12:00:00", scale="utc")

# Set up the plot
fig, ax = plt.subplots(figsize=(10, 10))

# Add interactive tools
fig.canvas.manager.toolbar.pan()
fig.canvas.manager.toolbar.zoom()

# Enable mousewheel zooming
def zoom_func(event):
    base_scale = 1.5
    if event.button == 'down':
        scale_factor = 1 / base_scale
    elif event.button == 'up':
        scale_factor = base_scale
    else:
        return

    x, y = event.xdata, event.ydata
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    xdata_range = xlim[1] - xlim[0]
    ydata_range = ylim[1] - ylim[0]
    xdata = xlim[0] + xdata_range * (x - xlim[0]) / xdata_range
    ydata = ylim[0] + ydata_range * (y - ylim[0]) / ydata_range
    ax.set_xlim([x - xdata_range / (2 * scale_factor), x + xdata_range / (2 * scale_factor)])
    ax.set_ylim([y - ydata_range / (2 * scale_factor), y + ydata_range / (2 * scale_factor)])
    fig.canvas.draw_idle()

fig.canvas.mpl_connect('scroll_event', zoom_func)

# Add the animation button
ax_button = plt.axes([0.81, 0.01, 0.1, 0.05])
button = Button(ax_button, 'Animate')
button.on_clicked(animate)

# Enable interactive mode
plt.ion()

# Plot the initial solar system
plot_solar_system(t)

# Show the plot
plt.show()

# Keep the plot open until it is manually closed
while True:
    plt.pause(0.1)
    if not plt.get_fignums():
        breaks