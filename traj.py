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

    r = np.random.random(2)

    x0 = np.array([r[0], 1-r[0]])
    y0 = np.array([r[1], 1-r[1]])

    T_vals = (0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.5, 2.0)
    fig, axs = plt.subplots(3, 3)
    fig.tight_layout()
    
    for i, T in enumerate(T_vals):
        ax = axs[i//3, i%3]
    
        sys = ql.QLD(A, A.T, T, T)
        x, y, t = sys.simulate(0.1, x0=1.0*x0, y0=1.0*y0)
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
        ax.set_title(rf"$T = {T}$, converges to $({signif(x_p[-1], 2)}, {signif(y_p[-1], 2)})$", fontsize=10)
        cb.ax.set_title(r"$t$")
        ax.scatter(r[0], r[1], color='#440256')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel(r"$x_1$")
        ax.set_ylabel(r"$y_1$")
        ax.set_aspect('equal')


    plt.show()


q6_1()