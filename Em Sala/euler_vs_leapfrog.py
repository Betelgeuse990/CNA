"""
Cálculo Numérico Aplicado - Aula de 02/07/2026
Método de Euler vs Leapfrog para EDOs

Contexto: simulação de interação gravitacional entre dois corpos.
"""
import numpy as np
import matplotlib.pyplot as plt

# Inputs
G = 1.0     # gravitational constant

m1 = 1.0    # [kg] mass of body 1
m2 = 1.0    # [kg] mass of body 2

r1 = np.array([-0.5, 0.0])    # position vectors
r2 = np.array([0.5, 0.0])

v1 = np.array([0.0, -0.5])    # initial velocities
v2 = np.array([0.0, 0.5])

t0 = 0.0        # simulation time start
tf = 100.0      # sim. end
delta_t = 1e-3  # time step

n = int((tf - t0) / delta_t)    # number of points in the grid

# Acceleration
def acc(r_1: np.ndarray, m_1, r_2: np.ndarray, m_2) -> tuple:
    """
    Calculates the acceleration vectors for Force and Body Mass
    of both bodies. Returns (vec_acc_b1, vec_acc_b2)

    d²r/dt² = F / m
    """
    dr = r_2 - r_1
    dist = np.linalg.norm(dr)

    a1 = G * m_2 * dr / dist ** 3
    a2 = -G * m_1 * dr / dist ** 3
    
    return a1, a2


def integrate_euler(r1_0, r2_0, v1_0, v2_0, t0, tf, delta_t) -> tuple:
    """
        dr/dt = v
        dv/dt = a
        
        Euler:
            r[i+1] = r[i] + delta_t * v[i]
            v[i+1] = v[i] + delta_t * a[i]
    """
    n = int((tf - t0) / delta_t)

    # Time array
    ts = np.zeros(n + 1)

    # Position history
    r1 = np.zeros((n + 1, 2))
    r2 = np.zeros((n + 1, 2))

    # Velocity history
    v1 = np.zeros((n + 1, 2))
    v2 = np.zeros((n + 1, 2))

    # Initial conditions
    r1[0] = r1_0
    r2[0] = r2_0

    v1[0] = v1_0
    v2[0] = v2_0

    for i in range(n):
        ts[i + 1] = ts[i] + delta_t

        a1, a2 = acc(r1[i], m1, r2[i], m2)

        r1[i + 1] = r1[i] + delta_t * v1[i]
        r2[i + 1] = r2[i] + delta_t * v2[i]

        v1[i + 1] = v1[i] + delta_t * a1
        v2[i + 1] = v2[i] + delta_t * a2

    return ts, r1, r2, v1, v2
    

def integrate_leapfrog(r1_0, r2_0, v1_0, v2_0, t0, tf, delta_t) -> tuple:
    """
    Leapfrog / velocity-Verlet integrator for the gravitational two-body problem.

    Updates:
        v_half = v[i] + 0.5 * dt * a[i]
        r[i+1] = r[i] + dt * v_half
        v[i+1] = v_half + 0.5 * dt * a[i+1]

    In other terms: it gives a "half step" towards what Euler would use as v[i+1]
    but v[i+1] here uses this half step as its origin and a new acceleration based
    on these half-step position vectors.
    """
    n = int((tf - t0) / delta_t)

    # Time array
    ts = np.zeros(n + 1)
    ts[0] = t0

    # Position history
    r1 = np.zeros((n + 1, 2))
    r2 = np.zeros((n + 1, 2))

    # Velocity history
    v1 = np.zeros((n + 1, 2))
    v2 = np.zeros((n + 1, 2))

    # Initial conditions
    r1[0] = r1_0
    r2[0] = r2_0

    v1[0] = v1_0
    v2[0] = v2_0

    for i in range(n):
        ts[i + 1] = ts[i] + delta_t

        a1, a2 = acc(r1[i], m1, r2[i], m2)

        # The Leapfrog half step
        v1_half = v1[i] + 0.5 * delta_t * a1
        v2_half = v2[i] + 0.5 * delta_t * a2

        r1[i + 1] = r1[i] + delta_t * v1_half
        r2[i + 1] = r2[i] + delta_t * v2_half

        # New acceleration vectors
        new_a1, new_a2 = acc(r1[i + 1], m1, r2[i + 1], m2)
        
        v1[i + 1] = v1_half + 0.5 * delta_t * new_a1
        v2[i + 1] = v2_half + 0.5 * delta_t * new_a2

    return ts, r1, r2, v1, v2


# ---------------------------------------------------------------------------------
def mechanical_energy(r1, r2, v1, v2):
    """
    Computes total mechanical energy E = K + U at each time step.
    """
    dist = np.linalg.norm(r2 - r1, axis=1)

    K = 0.5 * m1 * np.sum(v1**2, axis=1) + 0.5 * m2 * np.sum(v2**2, axis=1)
    U = -G * m1 * m2 / dist

    return K + U


def save_results(filename, ts, r1, r2, v1, v2, E):
    data = np.column_stack([
        ts,
        r1[:, 0], r1[:, 1],
        v1[:, 0], v1[:, 1],
        r2[:, 0], r2[:, 1],
        v2[:, 0], v2[:, 1],
        E
    ])

    header = "t x1 y1 vx1 vy1 x2 y2 vx2 vy2 E"

    np.savetxt(filename, data, header=header, comments="", fmt="%.10e")


def plot_trajectory(filename, r1, r2, title):
    plt.figure(figsize=(6, 6))

    plt.plot(r1[:, 0], r1[:, 1], label="Body 1")
    plt.plot(r2[:, 0], r2[:, 1], label="Body 2")

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(title)
    plt.axis("equal")
    plt.grid(True)
    plt.legend()

    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()


def plot_phase_space(filename, r1, r2, v1, v2, title):
    plt.figure(figsize=(8, 5))

    plt.plot(r1[:, 0], v1[:, 0], label="Body 1: x1 vs vx1")
    plt.plot(r2[:, 0], v2[:, 0], label="Body 2: x2 vs vx2")

    plt.xlabel("x")
    plt.ylabel("vx")
    plt.title(title)
    plt.grid(True)
    plt.legend()

    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()


def plot_energy_comparison(filename, ts_euler, E_euler, ts_leapfrog, E_leapfrog):
    plt.figure(figsize=(8, 5))

    plt.plot(ts_euler, E_euler, label="Euler")
    plt.plot(ts_leapfrog, E_leapfrog, label="Leapfrog")

    plt.xlabel("t")
    plt.ylabel("Total mechanical energy E")
    plt.title("Total energy: Euler vs Leapfrog")
    plt.grid(True)
    plt.legend()

    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()


# Run simulation
if __name__ == "__main__":
    # Euler
    ts_euler, r1_euler, r2_euler, v1_euler, v2_euler = integrate_euler(
        r1, r2, v1, v2, t0, tf, delta_t
    )

    E_euler = mechanical_energy(r1_euler, r2_euler, v1_euler, v2_euler)

    save_results(
        "euler.dat",
        ts_euler,
        r1_euler,
        r2_euler,
        v1_euler,
        v2_euler,
        E_euler
    )

    plot_trajectory(
        "euler_trajectory.png",
        r1_euler,
        r2_euler,
        "Euler trajectory"
    )

    plot_phase_space(
        "euler_phase_space.png",
        r1_euler,
        r2_euler,
        v1_euler,
        v2_euler,
        "Euler phase space"
    )

    # Leapfrog
    ts_leapfrog, r1_leapfrog, r2_leapfrog, v1_leapfrog, v2_leapfrog = integrate_leapfrog(
        r1, r2, v1, v2, t0, tf, delta_t
    )

    E_leapfrog = mechanical_energy(r1_leapfrog, r2_leapfrog, v1_leapfrog, v2_leapfrog)

    save_results(
        "leapfrog.dat",
        ts_leapfrog,
        r1_leapfrog,
        r2_leapfrog,
        v1_leapfrog,
        v2_leapfrog,
        E_leapfrog
    )

    plot_trajectory(
        "leapfrog_trajectory.png",
        r1_leapfrog,
        r2_leapfrog,
        "Leapfrog trajectory"
    )

    plot_phase_space(
        "leapfrog_phase_space.png",
        r1_leapfrog,
        r2_leapfrog,
        v1_leapfrog,
        v2_leapfrog,
        "Leapfrog phase space"
    )

    # Energy comparison
    plot_energy_comparison(
        "energy_comparison.png",
        ts_euler,
        E_euler,
        ts_leapfrog,
        E_leapfrog
    )

    print("Generated files:")
    print("  euler.dat")
    print("  leapfrog.dat")
    print("  euler_trajectory.png")
    print("  leapfrog_trajectory.png")
    print("  euler_phase_space.png")
    print("  leapfrog_phase_space.png")
    print("  energy_comparison.png")
    