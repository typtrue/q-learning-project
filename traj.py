import matplotlib.pyplot as plt 
import matplotlib.collections as col
import numpy as np
import ql_tools as ql

def signif(x, p):
    x = np.asarray(x)
    x_positive = np.where(np.isfinite(x) & (x != 0), np.abs(x), 10**(p-1))
    mags = 10 ** (p - 1 - np.floor(np.log10(x_positive)))
    return np.round(x * mags) / mags

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
        x_p = [k[0] for k in x]
        y_p = [k[0] for k in y]
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
    GRIDSIZE = 7
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
            xs[i] = x[-1][0]
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

def 

# q6_1()
# q6_2()