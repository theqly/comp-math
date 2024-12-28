import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

def read_data(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    times = []
    data = []
    current_time = None
    current_data = []

    for line in lines:
        if line.startswith("# Time ="):
            if current_time is not None:
                times.append(current_time)
                data.append(np.array(current_data))
            current_time = float(line.split("=")[1].strip())
            current_data = []
        elif line.startswith("#"):
            continue
        elif line.strip():
            current_data.append(list(map(float, line.strip().split())))

    if current_time is not None:
        times.append(current_time)
        data.append(np.array(current_data))

    return times, data

def animate_results(filename, save_as_gif=False):
    times, data = read_data(filename)

    x = data[0][:, 0]
    exact = [frame[:, 1] for frame in data]
    godunov = [frame[:, 2] for frame in data]
    implicit = [frame[:, 3] for frame in data]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlabel("x")
    ax.set_ylabel("u(t, x)")
    ax.set_title("Evolution of Numerical and Exact Solutions")

    exact_line, = ax.plot([], [], 'k-', label='Exact solution')
    godunov_line, = ax.plot([], [], 'ro-', label='Godunov scheme')
    implicit_line, = ax.plot([], [], 'bs-', label='Implicit scheme')
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

    ax.legend()

    def init():
        exact_line.set_data([], [])
        godunov_line.set_data([], [])
        implicit_line.set_data([], [])
        time_text.set_text('')
        return exact_line, godunov_line, implicit_line, time_text

    def update(frame):
        exact_line.set_data(x, exact[frame])
        godunov_line.set_data(x, godunov[frame])
        implicit_line.set_data(x, implicit[frame])
        time_text.set_text(f'Time = {times[frame]:.2f}')
        return exact_line, godunov_line, implicit_line, time_text

    ani = FuncAnimation(fig, update, frames=len(times), init_func=init, blit=True, interval=200)

    if save_as_gif:
        gif_writer = PillowWriter(fps=5)
        ani.save("animation.gif", writer=gif_writer)
        print("GIF saved as animation.gif")
    else:
        plt.show()

if __name__ == "__main__":
    animate_results("results.dat", save_as_gif=True)
