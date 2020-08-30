from functools import partial

import numpy as np
import plotly.graph_objects as go
from sympy import symbols
from sympy.utilities.lambdify import lambdify

from space_sim.maths.constants import G
from space_sim.maths.workflow import main


def _make_scatter(fn, axis, base_name, x_vals, line_shape):
    return go.Scatter(
        x=list(x_vals),
        y=[fn(x_val) for x_val in x_vals],
        line_shape=line_shape,
        name=f"{base_name} - {axis}",
    )


def make_plot(T, fn, title, base_name, line_shape):
    t_vals = np.linspace(0, float(T), 99)

    scatter = partial(
        _make_scatter, base_name=base_name, x_vals=t_vals, line_shape=line_shape
    )

    fig = go.Figure()
    fig.add_trace(scatter(partial(fn, i=0), "x"))
    fig.add_trace(scatter(partial(fn, i=1), "y"))
    fig.add_trace(scatter(partial(fn, i=2), "z"))

    fig.update_layout(title=title, xaxis_title="t", yaxis_title=base_name)
    return fig


def evaluate(fn, symbols, values):
    association = {symbols[key]: values[key] for key in symbols.keys()}
    return fn.subs({**association, G: 9.807})


def to_lambda(t, fn):
    lambdified = lambdify(t, fn)

    def final_fn(t, i):
        return lambdified(t)[i][0]

    return final_fn


def plot_position(t, T, p_fn, frames, variables, **args):
    T = evaluate(T, variables, args)
    p_fn = to_lambda(t, evaluate(p_fn.to_matrix(frames["N"]), variables, args))

    plot = make_plot(T, p_fn, "Position", "p(t)", line_shape="spline")
    plot.show()


def plot_velocity(t, T, v_fn, frames, variables, **args):
    T = evaluate(T, variables, args)
    v_fn = to_lambda(t, evaluate(v_fn.to_matrix(frames["N"]), variables, args))

    plot = make_plot(T, v_fn, "Velocity", "v(t)", line_shape="spline")
    plot.show()


def plot_acceleration(t, T, a_fn, frames, variables, **args):
    T = evaluate(T, variables, args)
    a_fn = to_lambda(t, evaluate(a_fn.to_matrix(frames["N"]), variables, args))

    plot = make_plot(T, a_fn, "Acceleration", "a(t)", line_shape="hv")
    plot.show()


def make_scatter3d(base_name, t_vals, p, v, a, p_mag, v_mag, a_mag):
    scatter3d = go.Scatter3d(
        x=[p(t_val, 0) for t_val in t_vals],
        y=[p(t_val, 1) for t_val in t_vals],
        z=[p(t_val, 2) for t_val in t_vals],
        line={
            "color": list(t_vals),
            "colorscale": "Viridis",
            "colorbar": {"title": "t"},
        },
        marker={"size": 2, "color": list(t_vals), "colorscale": "Viridis"},
        name=f"{base_name}(t)",
    )
    scatter3d.customdata = np.dstack(
        (
            [p_mag(t_val) for t_val in t_vals],
            [v_mag(t_val) for t_val in t_vals],
            [v(t_val, 0) for t_val in t_vals],
            [v(t_val, 1) for t_val in t_vals],
            [v(t_val, 2) for t_val in t_vals],
            [a_mag(t_val) for t_val in t_vals],
            [a(t_val, 0) for t_val in t_vals],
            [a(t_val, 2) for t_val in t_vals],
            [a(t_val, 2) for t_val in t_vals],
        )
    )[0]
    scatter3d.hovertemplate = """
    p:%{customdata[0]:.1f} m - (x:%{x:.1f}, y:%{y:.1f}, z:%{z:.1f})<br>
    v:%{customdata[1]:.1f} m/s - (x:%{customdata[2]:.1f}, y:%{customdata[3]:.1f}, z:%{customdata[4]:.1f})<br>
    a:%{customdata[5]:.1f} m/s2 - (x:%{customdata[6]:.1f}, y:%{customdata[7]:.1f}, z:%{customdata[8]:.1f})
    """
    return scatter3d


def make_3d_plot(T, p, v, a, p_mag, v_mag, a_mag, title, base_name):
    t_vals = np.linspace(0, float(T), 99)

    fig = go.Figure()
    fig.add_trace(make_scatter3d(base_name, t_vals, p, v, a, p_mag, v_mag, a_mag))

    fig.update_layout(
        title=title,
        xaxis_title=f"{base_name}(x)",
        yaxis_title=f"{base_name}(y)",
        autosize=True,
        height=700,
    )
    return fig


def plot_3d_position(t, T, a_fn, v_fn, p_fn, frames, variables, **args):
    T = evaluate(T, variables, args)
    a_fn_N = to_lambda(t, evaluate(a_fn.to_matrix(frames["N"]), variables, args))
    v_fn_N = to_lambda(t, evaluate(v_fn.to_matrix(frames["N"]), variables, args))
    p_fn_N = to_lambda(t, evaluate(p_fn.to_matrix(frames["N"]), variables, args))

    a_mag = lambdify(t, evaluate(a_fn.magnitude(), variables, args))
    v_mag = lambdify(t, evaluate(v_fn.magnitude(), variables, args))
    p_mag = lambdify(t, evaluate(p_fn.magnitude(), variables, args))

    plot = make_3d_plot(
        T, p_fn_N, v_fn_N, a_fn_N, p_mag, v_mag, a_mag, "3D Position", "p"
    )
    plot.show(config={"scrollZoom": False})


def wrap_plots_in_system(*plots):
    t = symbols("t", real=True, nonnegative=True)
    T, a_fn, v_fn, p_fn, frames, variables = main(t)

    return [
        partial(
            plot,
            t=t,
            T=T,
            a_fn=a_fn,
            v_fn=v_fn,
            p_fn=p_fn,
            frames=frames,
            variables=variables,
        )
        for plot in plots
    ]
