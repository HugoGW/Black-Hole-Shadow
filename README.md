# Photon Trajectories around a Schwarzschild Black Hole

This project simulates and animates the bending of light (null geodesics) around a non-rotating (Schwarzschild) black hole. It visualizes **gravitational lensing** and the formation of the **black hole shadow** using a wall of photons launched toward the black hole. The simulation includes:

- The **event horizon** (black disk)
- The **photon sphere** (unstable circular orbit for photons)
- The **shadow radius** (apparent dark region caused by gravitational deflection)
- An animated stream of photons interacting with the gravitational field

## What You See

- **Yellow lines** represent individual photon trajectories.
- Photons passing close to the black hole bend significantly.
- Photons with small enough impact parameters spiral into the black hole.
- The **shadow** of the black hole appears as a dark circular region, bounded by the **critical impact parameter** leading to unstable orbits at the **photon sphere**.

---

## Mathematical Background

This simulation uses the geodesic equations derived from the Schwarzschild metric in general relativity for **null (photon-like) trajectories** in the equatorial plane.

### Schwarzschild Metric (in natural units):
$$ \displaystyle ds^2 = -\big(1 - \frac{R_s}{r} \big) dt^2 + \big(1 - \frac{R_s}{r} \big)^{-1} dr^2 + r^2 d\phi^2$$

For photons $( ds^2 = 0 ) $, the motion reduces to a system of ODEs describing radial and angular evolution over coordinate time $t$.

### Geodesic Equations (Equatorial Plane):
The code evolves the following second-order differential equations:

$$\displaystyle \frac{d^2 r}{dt^2} = r \big( \frac{d\phi}{dt} \big)^2 - \frac{3}{2} R_s \big( \frac{d\phi}{dt} \big)^2$$

$$\displaystyle \frac{d^2 \phi}{dt^2} = -\frac{2}{r} \frac{dr}{dt} \frac{d\phi}{dt}$$

These are solved using a **custom Runge-Kutta 4 (RK4)** integrator accelerated with **Numba** for performance.

---

## Physical Features

- **Photon Sphere**: Located at $r = 1.5 R_s$, this is where photons can orbit the black hole on unstable circular paths.
- **Black Hole Shadow**: The apparent radius of the shadow observed from infinity is:
- 
  $$R_{\text{shadow}} = \frac{\sqrt{27}}{2} R_s \approx 2.6 R_s$$
  
  This is the boundary between light that escapes to infinity and light that is captured.

---

## How It Works

1. A wall of photons is initialized with varying vertical offsets, all aimed toward the black hole from a fixed horizontal distance.
2. Each photon’s trajectory is integrated using the Schwarzschild geodesic equations.
3. If a photon crosses the **event horizon** $( r < R_s )$, its trajectory is terminated and faded out in the animation.
4. The photons’ deflections reveal the characteristic circular **shadow** of the black hole.
5. The **photon sphere** and **shadow boundary** are drawn for reference.

---

## Dependencies

- `numpy`
- `matplotlib`
- `numba`
- `scipy` (for initial testing or extension, not used in the final animation)
- `ffmpeg` (optional, for saving animations as MP4)

You can install the required Python packages with:

```bash
pip install numpy matplotlib numba
````

To save the animation:

```python
anim.save("black_hole_shadow.mp4", writer='ffmpeg', fps=30, dpi=200)
```

---

## Purpose and Application

This simulation is designed to **visualize the gravitational effects** of a black hole on light rays. It illustrates how the **black hole shadow** forms and allows exploration of:

* Gravitational lensing
* The geometry of spacetime around compact objects
* The observable size of black holes (important for EHT imagery like M87\*)

---

## Customization

You can control:

* `Rs`: Schwarzschild radius
* `n_photons`: Number of photons in the wall
* `x0`: Initial distance from the black hole
* `fade_duration`: Fade-out duration of absorbed photons
* Appearance: Colors, background, and overlay visuals (e.g., photon sphere, shadow)

---

## Possible Extensions

* Add rotating (Kerr) black hole support
* Trace massive particles (timelike geodesics)
* Simulate gravitational lensing with background stars
* Include redshift calculation for each photon path

---

## References

* **Misner, Thorne, Wheeler** – *Gravitation* (1973)
* **Carroll, Sean** – *Spacetime and Geometry* (2004)
* Event Horizon Telescope collaboration papers
* Numerical methods: RK4 integration, Numba acceleration

---

Enjoy exploring the warped paths of photons in curved spacetime!
