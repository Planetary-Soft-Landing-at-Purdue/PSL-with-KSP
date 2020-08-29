from scipy.optimize import newton
from numpy import sin as s, cos as c, log, power as pow, sqrt, pi
import numpy as np
import krpc

def greatCircleDistance(coord1, coord2, r):
    Φ1 = np.deg2rad(coord1[0]); Φ2 = np.deg2rad(coord2[0])
    λ1 = np.deg2rad(coord1[1]); λ2 = np.deg2rad(coord2[1])
    Δσ = np.arccos(s(Φ1)*s(Φ2) + c(Φ1)*c(Φ2)*c(abs(λ2-λ1)))
    return r*Δσ

def E(t): return log(1-md*t/m0)
def F(t): return (1-md*t/m0)*(E(t)-1)+1

def x(θ, t): return x0 + u0*t + m0*T*c(θ)*F(t)/pow(md,2)
def h(θ, t): return h0 + v0*t - 0.5*g*pow(t,2) + m0*T*s(θ)*F(t)/pow(md,2)
def u(θ, t): return u0 - T*c(θ)*E(t)/md
def v(θ, t): return v0 - g*t - T*s(θ)*E(t)/md

def hθ(θ, t): return m0*T*c(θ)*F(t)/pow(md,2)
def uθ(θ, t): return T*s(θ)*E(t)/md
def vθ(θ, t): return -T*c(θ)*E(t)/md

def func_θ(θc, tf):
    hf = h(θc, tf)
    uf = u(θc, tf)
    vf = v(θc, tf)
    return (
        uθ(θc,tf) * (vf*sqrt(pow(vf,2)+2*g*hf) + pow(vf,2) + 2*g*hf) +
        uf*(vθ(θc,tf)*(vf+sqrt(pow(vf,2)+2*g*hf)) + g*hθ(θc,tf))
    )

def x_max(θc, tf):
    hf = h(θc, tf)
    uf = u(θc, tf)
    vf = v(θc, tf)
    h_pdg = 0
    return -(uf/g) * (vf + sqrt(pow(vf,2) - 2*g*(h_pdg-hf)))


# update state
def update_state(vessel, pcpf):
    body = vessel.orbit.body
    coord = vessel.flight().latitude, vessel.flight().longitude
    coord_pad = (body.latitude_at_position((0,0,0), pcpf),
                 body.longitude_at_position((0,0,0), pcpf))
    x0 = greatCircleDistance(coord, coord_pad, body.equatorial_radius)
    h0 = vessel.flight().mean_altitude
    u0 = vessel.flight().horizontal_speed(body.reference_frame)
    v0 = vessel.flight().vertical_speed(body.reference_frame)
    m0 = vessel.mass
    return x0, h0, u0, v0, m0

# secant search optimal thrust angle and t-cutoff
def optimize_bbg(tf1, tf2):
    θc1 = newton(func_θ, pi, args=(tf1,))
    θc2 = newton(func_θ, pi, args=(tf2,))
    xm1 = x_max(θc1, tf1)
    xm2 = x_max(θc2, tf2)
    while abs(tf2 - tf1) > 0.001:
        tf3 = tf2 - ((tf2 - tf1) / (xm2 - xm1)) * (xm2 - x0)
        tf1 = tf2; tf2 = tf3
        θc1 = newton(func_θ, pi, args=(tf1,))
        θc2 = newton(func_θ, pi, args=(tf2,))
        xm1 = x_max(θc1, tf1)
        xm2 = x_max(θc2, tf2)
    return θc2, tf2


# makes all necessary connections with ksp
conn = krpc.connect(address='192.168.0.109')
sc = conn.space_center
vessel = sc.active_vessel
ap = vessel.auto_pilot
orbit = vessel.orbit
body = orbit.body
bcbf = body.reference_frame
g = body.surface_gravity

T = vessel.max_vacuum_thrust
isp = vessel.vacuum_specific_impulse
md = T/(isp*9.81)
tf1 = 1; tf2 = 2
initial_heading = vessel.flight().heading
ap.engage()

while tf > 1:
    x0, h0, u0, v0, m0 = update_state()
    θc, tf = optimize_bbg(tf1, tf2)
    ap.target_heading = (initial_heading + 180) % 360
    ap.target_pitch = np.rad2deg(θc) - 180
    ap.target_roll = 0.0
    vessel.control.throttle = 1.0

met_1 = vessel.met
while vessel.met < met_1 + 1:
    pass

vessel.control.throttle = 0.0