import numpy as np
import scipy.integrate as sci

class QLD(object):
    """j"""
    def __init__(self, A: np.ndarray, B: np.ndarray, T_x : float, T_y : float, *, verbose=True) -> None:
        self.mats = (A, B)
        self.T = (T_x, T_y)
        self.verbose = verbose

    def qld_eq(self, player : int, x : np.ndarray, y : np.ndarray) -> np.ndarray:
        """
        Return time-derivative of player's choices using the Q-Learning Dynamic.
        
        Parameters
        ----------
        player : int
            Player to calculate time derivative for. Must be 0 or 1.
        x : np.ndarray
            Player X's strategy at current timestep.
        y : np.ndarray
            Player Y's strategy at current timestep.
        
        Returns
        -------
        res : np.ndarray
            Time derivative of player's choices at current timestep.
        """

        s = (x, y)
        mats = self.mats
        T = self.T

        if not np.isclose(sum(x), 1):
            raise ValueError(f"x must be a probability vector, got {x}")
        if not np.isclose(sum(y), 1):
            raise ValueError(f"y must be a probability vector, got {y}")
        if (player != 0 and player != 1):
            raise ValueError(f"player must be 0 or 1, got {player}")
        if len(x) != mats[0].shape[1]:
            raise ValueError(f"""length of x and number of columns of A must 
            match, got {len(x)} and {mats[0].shape[0]} respectively""")
        if len(y) != mats[0].shape[0]:
            raise ValueError(f"""length of y and number of rows of A must 
            match, got {len(y)} and {mats[0].shape[1]} respectively""")

        n = len(s[player])

        # p1 = Ay or x^T B depending on player
        if player:
            p1 = s[1-player].T @ mats[player]
        else:
            p1 = mats[player] @ s[1-player]

        # p2_i = -x^T M_p y + T_p sum_{j=1}^n s_j ln(s_j / s_i)
        p2 = np.array([-x.T @ mats[player] @ y + T[player] * sum([sj * np.log(sj / si) for sj in s[player]]) for si in s[player]])

        if any(np.isnan(p2)):
            print(f"s: {s[player]}")
            raise ValueError

        return s[player] * (p1 + p2)
    
    def simulate(self, dt : float, *, maxits : int = 10000, tol : float = 10**(-6), x0 : np.ndarray = None, y0 : np.ndarray = None):
        """
        
        """
        n = self.mats[0].shape[0]

        if x0 is None:
            x0 = np.ones(n[1]) / n[1]
        if y0 is None:
            y0 = np.ones(n[0]) / n[0]

        xs = np.zeros((n, maxits+1))
        ys = np.zeros((n, maxits+1))

        xs[:, 0] = x0
        ys[:, 0] = y0
        ts = [0]
        t = 0

        for i in range(maxits):
            dx = self.qld_eq(0, x0, y0)
            dy = self.qld_eq(1, x0, y0)

            if np.linalg.norm(dx) < tol and np.linalg.norm(dy) < tol:
                if self.verbose:
                    print(f"converged after {i} iterations!")
                break

            x0 += dx * dt
            y0 += dy * dt

            xs[:, i+1] = x0
            ys[:, i+1] = y0

            t += dt
            ts.append(t)
            if (i + 1) % 100 == 0 and self.verbose:
                print(f"{i+1}/{maxits}")

        return xs[:, :i+1], ys[:, :i+1], ts

    def equilibria(self):
        raise NotImplementedError