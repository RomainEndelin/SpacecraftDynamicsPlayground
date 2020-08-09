def v_args(args):
    return {
        a_magnitude: args["a_knob"],
        v0_x: args["v0x_knob"],
        v0_y: args["v0y_knob"],
        v0_z: args["v0z_knob"],
        vT_x: args["vTx_knob"],
        vT_y: args["vTy_knob"],
        vT_z: args["vTz_knob"],
    }

def p_args(args):
    return {
        **v_args(args),
        p0_x: args["p0x_knob"],
        p0_y: args["p0y_knob"],
        p0_z: args["p0z_knob"],
        pT_x: args["pTx_knob"],
        pT_y: args["pTy_knob"],
        pT_z: args["pTz_knob"],
    }

def make_plot(title, short_name, fn, args, line_shape="spline"):
    R_val = solved_R.subs(v_args(args))
    T_val = solved_T.subs(p_args(args))

    t_sample = np.linspace(0, float(T_val), 99)
    image_sample = lambdify(t, fn.to_matrix(N).subs({**p_args(args), T: T_val, R: R_val}))

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=list(t_sample), y=[image_sample(t_val)[0][0] for t_val in t_sample], line_shape=line_shape, name=f"{short_name} - x"))
    fig.add_trace(go.Scatter(x=list(t_sample), y=[image_sample(t_val)[1][0] for t_val in t_sample], line_shape=line_shape, name=f"{short_name} - y"))
    fig.add_trace(go.Scatter(x=list(t_sample), y=[image_sample(t_val)[2][0] for t_val in t_sample], line_shape=line_shape, name=f"{short_name} - z"))

    fig.update_layout(title=title, xaxis_title="t", yaxis_title=short_name)
    return fig

def plot_position(**knob_args):
    plot = make_plot("Position", "p(t)", d_fn, knob_args)
    plot.show()

def plot_velocity(**knob_args):
    plot = make_plot("Velocity", "v(t)", v_fn, knob_args)
    plot.show()

def plot_acceleration(**knob_args):
    plot = make_plot("Acceleration", "a(t)", a_fn, knob_args, line_shape="hv")
    plot.show()

def plot_3d_position(**args):
    R_val = solved_R.subs(v_args(args))
    T_val = solved_T.subs(p_args(args))

    t_sample = np.linspace(0, float(T_val), 99)
    d = lambdify(t, d_fn.to_matrix(N).subs({**p_args(args), T: T_val, R: R_val}))
    p_mag = lambdify(t, d_fn.to_matrix(A).subs({**p_args(args), T: T_val, R: R_val}))
    v = lambdify(t, v_fn.to_matrix(N).subs({**p_args(args), T: T_val, R: R_val}))
    v_mag = lambdify(t, v_fn.to_matrix(A).subs({**p_args(args), T: T_val, R: R_val}))
    a = lambdify(t, a_fn.to_matrix(N).subs({**p_args(args), T: T_val, R: R_val}))
    a_mag = lambdify(t, a_fn.to_matrix(A).subs({**p_args(args), T: T_val, R: R_val}))

    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
        x=[d(t_val)[0][0] for t_val in t_sample],
        y=[d(t_val)[1][0] for t_val in t_sample],
        z=[d(t_val)[2][0] for t_val in t_sample],
        customdata=np.dstack((
            [p_mag(t_val)[0][0] for t_val in t_sample],
            [v_mag(t_val)[0][0] for t_val in t_sample],
            [v(t_val)[0][0] for t_val in t_sample],
            [v(t_val)[1][0] for t_val in t_sample],
            [v(t_val)[2][0] for t_val in t_sample],
            [a_mag(t_val)[0][0] for t_val in t_sample],
            [a(t_val)[0][0] for t_val in t_sample],
            [a(t_val)[1][0] for t_val in t_sample],
            [a(t_val)[2][0] for t_val in t_sample],
        ))[0],
        hovertemplate='''
            p:%{customdata[0]:.1f} m - (%{x:.1f} %{y:.1f} %{z:.1f})<br>
            v:%{customdata[1]:.1f} m/s - (x:%{customdata[2]:.1f}, y:%{customdata[3]:.1f}, z:%{customdata[4]:.1f})<br>
            a:%{customdata[5]:.1f} m/s2 - (x:%{customdata[6]:.1f}, y:%{customdata[7]:.1f}, z:%{customdata[8]:.1f})
        ''',
        line={
            "color": list(t_sample),
            "colorscale": "Viridis",
            "colorbar": { "title": "t" }
        },
        marker={ "size": 2, "color": list(t_sample), "colorscale": "Viridis" },
        name=f"p(t)",

    ))

    fig.update_layout(title="3D Position", xaxis_title="p(x)", yaxis_title="p(y)", autosize=True, height=700)
    fig.show(config={'scrollZoom': False})