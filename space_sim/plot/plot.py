from functools import partial

import numpy as np
import plotly.graph_objects as go
from sympy.utilities.lambdify import lambdify

from space_sim.maths.helpers import G
from space_sim.maths.main import (
    A,
    N,
    R,
    T,
    a_fn,
    a_magnitude,
    distance,
    p0_vector_N,
    p0x_n,
    p0y_n,
    p0z_n,
    p_fn,
    pT_vector_N,
    pTx_n,
    pTy_n,
    pTz_n,
    solved_R,
    solved_T,
    t,
    v0_magnitude,
    v_fn,
    vT_magnitude,
)


def v_args(args):
    return {
        G: 9.807,
        a_magnitude: args["a_magnitude"],
        v0_magnitude: args["v0_magnitude"],
        vT_magnitude: args["vT_magnitude"],
    }


def p_args(args):
    return {
        **v_args(args),
        p0x_n: args["p0x_n"],
        p0y_n: args["p0y_n"],
        p0z_n: args["p0z_n"],
        pTx_n: args["pTx_n"],
        pTy_n: args["pTy_n"],
        pTz_n: args["pTz_n"],
    }


def _make_scatter(fn, axis, base_name, x_vals, line_shape):
    return go.Scatter(
        x=list(x_vals),
        y=[fn(x_val) for x_val in x_vals],
        line_shape=line_shape,
        name=f"{base_name} - {axis}",
    )


def make_plot(title, base_name, fn, args, line_shape="spline"):
    d = (pT_vector_N - p0_vector_N).to_matrix(A).subs(p_args(args))[0]

    R_val = solved_R.subs(v_args(args))
    T_val = solved_T.subs({**p_args(args), distance: d})

    t_vals = np.linspace(0, float(T_val), 99)
    final_fn = lambdify(t, fn.to_matrix(N).subs({**p_args(args), T: T_val, R: R_val}),)

    scatter = partial(
        _make_scatter, base_name=base_name, x_vals=t_vals, line_shape=line_shape
    )

    fig = go.Figure()
    fig.add_trace(scatter(lambda t: final_fn(t)[0][0], "x"))
    fig.add_trace(scatter(lambda t: final_fn(t)[1][0], "y"))
    fig.add_trace(scatter(lambda t: final_fn(t)[2][0], "z"))

    fig.update_layout(title=title, xaxis_title="t", yaxis_title=base_name)
    return fig


def plot_position(**args):
    plot = make_plot("Position", "p(t)", p_fn + p0_vector_N, args)
    plot.show()


def plot_velocity(**args):
    plot = make_plot("Velocity", "v(t)", v_fn, args)
    plot.show()


def plot_acceleration(**args):
    plot = make_plot("Acceleration", "a(t)", a_fn, args, line_shape="hv")
    plot.show()


def plot_3d_position(**args):
    # TODO: Initial position
    d = (pT_vector_N - p0_vector_N).to_matrix(A).subs(p_args(args))[0]

    R_val = solved_R.subs(v_args(args))
    T_val = solved_T.subs({**p_args(args), distance: d})

    t_sample = np.linspace(0, float(T_val), 99)
    d = lambdify(t, p_fn.to_matrix(N).subs({**p_args(args), T: T_val, R: R_val}))
    p_mag = lambdify(t, p_fn.to_matrix(A).subs({**p_args(args), T: T_val, R: R_val}))
    v = lambdify(t, v_fn.to_matrix(N).subs({**p_args(args), T: T_val, R: R_val}))
    v_mag = lambdify(t, v_fn.to_matrix(A).subs({**p_args(args), T: T_val, R: R_val}))
    a = lambdify(t, a_fn.to_matrix(N).subs({**p_args(args), T: T_val, R: R_val}))
    a_mag = lambdify(t, a_fn.to_matrix(A).subs({**p_args(args), T: T_val, R: R_val}))

    fig = go.Figure()
    fig.add_trace(
        go.Scatter3d(
            x=[d(t_val)[0][0] for t_val in t_sample],
            y=[d(t_val)[1][0] for t_val in t_sample],
            z=[d(t_val)[2][0] for t_val in t_sample],
            customdata=np.dstack(
                (
                    [p_mag(t_val)[0][0] for t_val in t_sample],
                    [v_mag(t_val)[0][0] for t_val in t_sample],
                    [v(t_val)[0][0] for t_val in t_sample],
                    [v(t_val)[1][0] for t_val in t_sample],
                    [v(t_val)[2][0] for t_val in t_sample],
                    [a_mag(t_val)[0][0] for t_val in t_sample],
                    [a(t_val)[0][0] for t_val in t_sample],
                    [a(t_val)[1][0] for t_val in t_sample],
                    [a(t_val)[2][0] for t_val in t_sample],
                )
            )[0],
            hovertemplate="""
            p:%{customdata[0]:.1f} m - (x:%{x:.1f}, y:%{y:.1f}, z:%{z:.1f})<br>
            v:%{customdata[1]:.1f} m/s - (x:%{customdata[2]:.1f}, y:%{customdata[3]:.1f}, z:%{customdata[4]:.1f})<br>
            a:%{customdata[5]:.1f} m/s2 - (x:%{customdata[6]:.1f}, y:%{customdata[7]:.1f}, z:%{customdata[8]:.1f})
        """,
            line={
                "color": list(t_sample),
                "colorscale": "Viridis",
                "colorbar": {"title": "t"},
            },
            marker={"size": 2, "color": list(t_sample), "colorscale": "Viridis"},
            name="p(t)",
        )
    )

    fig.update_layout(
        title="3D Position",
        xaxis_title="p(x)",
        yaxis_title="p(y)",
        autosize=True,
        height=700,
    )
    fig.show(config={"scrollZoom": False})
