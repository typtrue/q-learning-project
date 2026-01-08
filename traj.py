import matplotlib.pyplot as plt 
import matplotlib.collections as col
import numpy as np
import ql_tools as ql

def signif(x, p):
    x = np.asarray(x)
    x_positive = np.where(np.isfinite(x) & (x != 0), np.abs(x), 10**(p-1))
    mags = 10 ** (p - 1 - np.floor(np.log10(x_positive)))
    return np.round(x * mags) / mags

def plot_helper(sys, x01, x02, x03, y01, y02, y03, filename, multistage=False):
    if multistage:
        # THIS IS A HORRIBLE WAY OF DOING THIS! FIX LATER
        xt1, yt1, _ = sys.simulate_multistage(0.1, maxits=200, x0=1.0*x01, y0=1.0*y01)
        xt2, yt2, _ = sys.simulate_multistage(0.1, maxits=200, x0=1.0*x02, y0=1.0*y02)
        xt3, yt3, _ = sys.simulate_multistage(0.1, maxits=200, x0=1.0*x03, y0=1.0*y03)
        plot_triangle(x01, y01, x03, xt1, yt1, xt3, filename, multistage)
    else:
        xt1, yt1, _ = sys.simulate(0.1, x0=1.0*x01, y0=1.0*y01)
        xt2, yt2, _ = sys.simulate(0.1, x0=1.0*x02, y0=1.0*y02)
        xt3, yt3, _ = sys.simulate(0.1, x0=1.0*x03, y0=1.0*y03)
        plot_triangle(x01, x02, x03, xt1, xt2, xt3, f"{filename}_01", multistage)
        plot_triangle(y01, y02, y03, yt1, yt2, yt3, f"{filename}_02", multistage)

def plot_triangle(x01, x02, x03, xt1, xt2, xt3, filename, multistage):
    # define a projection from the 3D simplex on a triangle
    proj = np.array(
    [[-1 * np.cos(30. / 360. * 2. * np.pi),np.cos(30. / 360. * 2. * np.pi),0.],
    [-1 * np.sin(30. / 360. * 2. * np.pi),-1 * np.sin(30. / 360. * 2. * np.pi),1.]])
    # project the boundary on the simplex onto the boundary of the triangle
    ts = np.linspace(0, 1, 10000)
    PBd1 = proj@np.array([ts,(1-ts),0*ts])
    PBd2 = proj@np.array([0*ts,ts,(1-ts)])
    PBd3 = proj@np.array([ts,0*ts,(1-ts)])

    # project the orbits on the triangle
    orbittriangle1=proj@xt1
    orbittriangle2=proj@xt2
    orbittriangle3=proj@xt3
    ic1=proj@x01
    ic2=proj@x02
    ic3=proj@x03
    # no box
    plt.box(False)
    plt.axis(False)
    # plot the orbits, the initial values, the corner points, and the boundary points
    plt.plot(orbittriangle1[0],orbittriangle1[1],"-",markersize=1,color='green')
    plt.plot(orbittriangle2[0],orbittriangle2[1],"-",markersize=1,color='red')
    plt.plot(ic1[0],ic1[1],"+",markersize=10,color='green')
    plt.plot(ic2[0],ic2[1],"+",markersize=10,color='red')
    if not multistage:
        plt.plot(orbittriangle3[0],orbittriangle3[1],"-",markersize=1,color='blue')
        plt.plot(ic3[0],ic3[1],"+",markersize=10,color='blue')
    plt.text(-0.8660254-0.2, -0.5 +0.05 , "rock",fontsize=12)
    plt.text(+0.8660254+0.05, -0.5 +0.05 , "paper",fontsize=12)
    plt.text(0-0.1, 1 +0.1 , "scissors",fontsize=12)
    plt.plot(PBd1[0], PBd1[1], ".",color='black',markersize=3)
    plt.plot(PBd2[0], PBd2[1], ".",color='black',markersize=3)
    plt.plot(PBd3[0], PBd3[1], ".",color='black',markersize=3)


    plt.savefig(f"{filename}.pdf")
    plt.show()

def q6_1():
    A = np.array([[6, 0], [4, 2]])

    X_SET = 0.483
    Y_SET = 0.820

    x0 = np.array([X_SET, 1-X_SET])
    y0 = np.array([Y_SET, 1-Y_SET])

    T_vals = (0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.5, 2.0)
    fig, axs = plt.subplots(3, 3)
    fig.tight_layout()
    fig.set_figheight(8)
    fig.set_figwidth(12)
    
    for i, T in enumerate(T_vals):
        ax = axs[i//3, i%3]
        ax.plot([0, 1], [0, 1], '--', c='black', linewidth=1, alpha=0.5)
    
        sys = ql.QLD(A, A.T, T, T)
        x, y, t = sys.simulate(0.1, x0=1.0*x0, y0=1.0*y0, maxits=100000)
        x_p = x[0, :]
        y_p = y[0, :]
        points = np.array([x_p, y_p]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        norm = plt.Normalize(min(t), max(t))
        lc = col.LineCollection(segments, cmap='viridis', norm=norm)
        lc.set_array(t)
        lc.set_linewidth(2)
        line = ax.add_collection(lc)
        cb = fig.colorbar(line, ax=ax)
        ax.set_title(rf"$T = {T}$, conv. to $({signif(x_p[-1], 2)}, {signif(y_p[-1], 2)})$", fontsize=10)
        cb.ax.set_title(r"$t$")
        ax.scatter(X_SET, Y_SET, color='#440256')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel(r"$x_1$")
        ax.set_ylabel(r"$y_1$")
        ax.set_aspect('equal')

    plt.savefig("subplots.pdf", format="pdf")
    plt.show()

def q6_2():
    RES = 50
    GRIDSIZE = 5
    EPS = 10**(-2)

    A = np.array([[6, 0], [4, 2]])

    r = np.random.random(2)

    fig, ax = plt.subplots()

    Ts = np.linspace(0, 1.5, RES)
    X = np.linspace(EPS, 1-EPS, GRIDSIZE)
    Y = np.linspace(EPS, 1-EPS, GRIDSIZE)
    for x_g, y_g in ((x, y) for x in X for y in Y):
        x0 = np.array([x_g, 1-x_g])
        y0 = np.array([y_g, 1-y_g])
        xs = np.zeros(RES)
        print(f"({x_g}, {y_g}): \n")
        for i, T in enumerate(Ts):
            sys = ql.QLD(A, A.T, T, T, verbose=False)
            x, _, t = sys.simulate(0.1, x0=1.0*x0, y0=1.0*y0)
            xs[i] = x[0, -1]
        points = np.array([Ts, xs]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        norm = plt.Normalize(min(Ts), max(Ts))
        lc = col.LineCollection(segments, cmap='viridis', norm=norm)
        lc.set_array(Ts)
        lc.set_linewidth(2)
        ax.add_collection(lc)
    ax.set_xlabel(r"$T$")
    ax.set_ylabel(r"$x_1$")
    ax.set_xlim(0, 1.5)
    ax.set_ylim(0, 1)
    plt.suptitle(r"""Bifurcation diagram of $x_1$ as $T$ increases""",
                    fontsize=18)
    ax.set_title(r"""Arrows drawn on to help illustrate stability""",
                 fontsize=8)

    plt.savefig("bifur.pdf", format="pdf")
    plt.show()


class RPS(ql.QLD):
    def __init__(self, eps_x: float, eps_y: float, b: np.ndarray | float, T_x : float, T_y : float, *, verbose=True) -> None:
        if isinstance(b, float):
            b_en = b * np.ones(3)
        else:
            b_en = 1.0 * b

        print(b_en)


        self.eps = (eps_x, eps_y)
        self.b = b_en

        A = np.array([[eps_x, -1, b_en[0]],
                      [b_en[1], eps_x, -1],
                      [-1, b_en[2], eps_x]])

        B = np.array([[eps_y, b_en[0], -1],
                      [-1, eps_y, b_en[1]],
                      [b_en[2], -1, eps_y]])

        super().__init__(A, B, T_x, T_y, verbose=verbose)

    def set_b(self, b: np.ndarray | float):
        if isinstance(b, float):
            b_en = b * np.ones(3)
        else:
            b_en = 1.0 * b
        eps_x, eps_y = self.eps

        A = np.array([[eps_x, -1, b_en[0]],
                      [b_en[1], eps_x, -1],
                      [-1, b_en[2], eps_x]])

        B = np.array([[eps_y, b_en[0], -1],
                      [-1, eps_y, b_en[1]],
                      [b_en[2], -1, eps_y]])

        self.mats = (A, B)

    def step(self, dt, x, y):
        dx = self.qld_eq(0, x, y)
        dy = self.qld_eq(1, x, y)

        x += dx * dt
        y += dy * dt

        return dx, dy

    def simulate_multistage(self, dt : float, *, maxits : int = 10000, tol : float = 10**(-6), x0 : np.ndarray = None, y0 : np.ndarray = None):
        """
        
        """
        n = self.mats[0].shape

        if x0 is None:
            x0 = np.ones(n[1]) / n[1]
        if y0 is None:
            y0 = np.ones(n[0]) / n[0]

        xs = np.zeros((3, maxits+1))
        ys = np.zeros((3, maxits+1))

        xs[:, 0] = x0
        ys[:, 0] = y0

        ts = [0]
        t = 0

        b = self.b

        for i in range(maxits):
            dx, dy = self.step(dt, x0, y0)

            if np.linalg.norm(dx) < tol and np.linalg.norm(dy) < tol and self.verbose:
                print(f"converged after {i} iterations!")
                break

            xs[:, i+1] = x0
            ys[:, i+1] = y0
            t += dt
            ts.append(t)

            rand_x = np.random.choice(3, p=x0)
            rand_y = np.random.choice(3, p=y0)
            
            if rand_x != rand_y:
                if rand_x % 3 == (rand_y + 1) % 3:
                    if b[rand_x] != 1:
                        b[rand_x] += 1
                    else:
                        b = np.ones(3)
                        b[rand_x] += 1
                elif rand_x % 3 == (rand_y - 1) % 3:
                    if b[rand_y] != 1:
                        b[rand_y] += 1
                    else:
                        b = np.ones(3)
                        b[rand_y] += 1
            
            self.set_b(b)

            if (i + 1) % 100 == 0 and self.verbose:
                print(f"{i+1}/{maxits}")

        return xs[:, :i+1], ys[:, :i+1], ts


def sec3():
    # todo: generate plots of how equilibria move from varying eps, b, and T respectively


    eps_def = 0.0
    b_def = 1.0
    T_def = 0.2

    x01 = np.array([0.01, 0.495, 0.495])
    x02 = np.array([0.495, 0.01, 0.495])
    x03 = np.array([0.495, 0.495, 0.01])

    y01 = np.array([0.495, 0.01, 0.495])
    y02 = np.array([0.495, 0.495, 0.01])
    y03 = np.array([0.01, 0.495, 0.495])

    sys = RPS(eps_def, eps_def, b_def, T_def, T_def)

    plot_helper(sys, x01, x02, x03, y01, y02, y03, "mod_rps")

def sec3_2():
    # todo: generate plots of how equilibria move from varying eps, b, and T respectively


    eps_def = 0.0
    b_def = 5.0
    T_def = 0.2

    x01 = np.array([0.01, 0.495, 0.495])
    x02 = np.array([0.495, 0.01, 0.495])
    x03 = np.array([0.495, 0.495, 0.01])

    y01 = np.array([0.495, 0.01, 0.495])
    y02 = np.array([0.495, 0.495, 0.01])
    y03 = np.array([0.01, 0.495, 0.495])

    sys = RPS(eps_def, eps_def, b_def, T_def, T_def)

    plot_helper(sys, x01, x02, x03, y01, y02, y03, "mod_rps_b_5")


def sec4():
    # todo: make game where successive wins with the same 

    eps_def = 0.0
    b_def = 1.0
    T_def = 0.5

    x01 = np.array([0.01, 0.495, 0.495])
    x02 = np.array([0.495, 0.01, 0.495])
    x03 = np.array([0.495, 0.495, 0.01])

    y01 = np.array([0.495, 0.01, 0.495])
    y02 = np.array([0.495, 0.495, 0.01])
    y03 = np.array([0.01, 0.495, 0.495])

    sys = RPS(eps_def, eps_def, b_def, T_def, T_def)

    plot_helper(sys, x01, x02, x03, y01, y02, y03, "multistage", multistage=True)

# q6_1()
# q6_2()
sec3()
#sec3_2()
# sec4()