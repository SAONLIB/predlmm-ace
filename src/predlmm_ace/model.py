"""Scalable PredLMM-ACE variance-component estimation.

The module never forms an N x N relationship, household, or inverse matrix.
It represents the genetic covariance as K = B B' + diag(delta), where B is
an on-disk or in-memory N x r Nystrom factor, and the shared-environment
covariance by integer family/household labels.

Two estimators are supplied:
  * ExactProfileML: exact profile maximum likelihood using Woodbury identities;
    appropriate when r is small enough that an r x r matrix can be factored.
  * MatrixFreeAIREML: matrix-free average-information REML with PCG solves and
    Hutchinson trace estimates; intended for very large cohorts such as UKB.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize
from scipy.sparse.linalg import LinearOperator, cg


FloatArray = NDArray[np.float64]


def _as_2d(x: ArrayLike, n: int) -> tuple[FloatArray, bool]:
    a = np.asarray(x, dtype=np.float64)
    was_vector = a.ndim == 1
    if was_vector:
        a = a[:, None]
    if a.ndim != 2 or a.shape[0] != n:
        raise ValueError(f"Expected ({n},) or ({n}, k), got {a.shape}.")
    return a, was_vector


def _restore_shape(x: FloatArray, was_vector: bool) -> FloatArray:
    return x[:, 0] if was_vector else x


class FamilyBlocks:
    """Implicit Pi matrix: Pi_ij=1 exactly when labels i and j match."""

    def __init__(self, labels: Sequence[object]):
        labels = np.asarray(labels)
        if labels.ndim != 1:
            raise ValueError("family_labels must be one-dimensional.")
        _, self.codes = np.unique(labels.astype(str), return_inverse=True)
        self.codes = self.codes.astype(np.int64, copy=False)
        self.n = self.codes.size
        self.n_groups = int(self.codes.max()) + 1 if self.n else 0
        self.sizes = np.bincount(self.codes, minlength=self.n_groups)

    def matmat(self, x: ArrayLike) -> FloatArray:
        a, was_vector = _as_2d(x, self.n)
        sums = np.zeros((self.n_groups, a.shape[1]), dtype=np.float64)
        np.add.at(sums, self.codes, a)
        return _restore_shape(sums[self.codes], was_vector)

    def base_inverse(
        self, x: ArrayLike, diagonal: ArrayLike, shared_variance: float
    ) -> FloatArray:
        """Apply (diag(diagonal) + shared_variance * Pi)^-1."""
        a, was_vector = _as_2d(x, self.n)
        diagonal = np.asarray(diagonal, dtype=np.float64)
        if diagonal.shape != (self.n,) or np.any(diagonal <= 0):
            raise ValueError("The base diagonal must be positive and length N.")
        if shared_variance < 0:
            raise ValueError("shared_variance must be nonnegative.")
        invd = 1.0 / diagonal
        y = invd[:, None] * a
        group_sum = np.zeros((self.n_groups, a.shape[1]), dtype=np.float64)
        np.add.at(group_sum, self.codes, y)
        group_invd = np.bincount(
            self.codes, weights=invd, minlength=self.n_groups
        )
        denom = 1.0 + shared_variance * group_invd
        correction = (
            shared_variance * invd[:, None] * group_sum[self.codes]
            / denom[self.codes, None]
        )
        return _restore_shape(y - correction, was_vector)

    def base_logdet(self, diagonal: ArrayLike, shared_variance: float) -> float:
        diagonal = np.asarray(diagonal, dtype=np.float64)
        if diagonal.shape != (self.n,) or np.any(diagonal <= 0):
            return np.inf
        invd = 1.0 / diagonal
        group_invd = np.bincount(
            self.codes, weights=invd, minlength=self.n_groups
        )
        return float(
            np.log(diagonal).sum()
            + np.log1p(shared_variance * group_invd).sum()
        )


class LowRankACE:
    """Implicit covariance components K_g, Pi, and I.

    B may be a numpy memmap. Matrix products are row-chunked so temporary
    storage is bounded independently of N.
    """

    def __init__(
        self,
        factor: ArrayLike,
        diagonal_correction: ArrayLike,
        family_labels: Sequence[object],
        row_chunk_size: int = 8192,
    ):
        self.B = np.asanyarray(factor)
        if self.B.ndim != 2:
            raise ValueError("factor must be an N x rank array.")
        self.n, self.rank = self.B.shape
        self.delta = np.asarray(diagonal_correction, dtype=np.float64)
        if self.delta.shape != (self.n,):
            raise ValueError("diagonal_correction must have length N.")
        if np.any(self.delta < -1e-8):
            raise ValueError("diagonal_correction cannot be materially negative.")
        self.delta = np.maximum(self.delta, 0.0)
        self.families = FamilyBlocks(family_labels)
        if self.families.n != self.n:
            raise ValueError("family_labels must have one label per sample.")
        self.row_chunk_size = int(row_chunk_size)

    def _rows(self):
        for start in range(0, self.n, self.row_chunk_size):
            yield slice(start, min(self.n, start + self.row_chunk_size))

    def bt_matmat(self, x: ArrayLike) -> FloatArray:
        a, was_vector = _as_2d(x, self.n)
        out = np.zeros((self.rank, a.shape[1]), dtype=np.float64)
        for sl in self._rows():
            out += np.asarray(self.B[sl], dtype=np.float64).T @ a[sl]
        return out[:, 0] if was_vector else out

    def b_matmat(self, x: ArrayLike) -> FloatArray:
        a = np.asarray(x, dtype=np.float64)
        was_vector = a.ndim == 1
        if was_vector:
            a = a[:, None]
        if a.ndim != 2 or a.shape[0] != self.rank:
            raise ValueError("The leading dimension must equal factor rank.")
        out = np.empty((self.n, a.shape[1]), dtype=np.float64)
        for sl in self._rows():
            out[sl] = np.asarray(self.B[sl], dtype=np.float64) @ a
        return _restore_shape(out, was_vector)

    def genetic(self, x: ArrayLike) -> FloatArray:
        a, was_vector = _as_2d(x, self.n)
        out = self.b_matmat(self.bt_matmat(a)) + self.delta[:, None] * a
        return _restore_shape(out, was_vector)

    def component(self, index: int, x: ArrayLike) -> FloatArray:
        if index == 0:
            return self.genetic(x)
        if index == 1:
            return self.families.matmat(x)
        if index == 2:
            return np.asarray(x, dtype=np.float64)
        raise IndexError("Component index must be 0 (genetic), 1 (family), or 2 (residual).")

    def covariance(self, theta: ArrayLike, x: ArrayLike) -> FloatArray:
        sg, sc, se = np.asarray(theta, dtype=np.float64)
        return sg * self.genetic(x) + sc * self.families.matmat(x) + se * np.asarray(x)

    def diagonal(self, theta: ArrayLike) -> FloatArray:
        sg, sc, se = np.asarray(theta, dtype=np.float64)
        row_norm = np.empty(self.n, dtype=np.float64)
        for sl in self._rows():
            b = np.asarray(self.B[sl], dtype=np.float64)
            row_norm[sl] = np.einsum("ij,ij->i", b, b)
        return sg * (row_norm + self.delta) + sc + se


@dataclass
class FitResult:
    genetic_variance: float
    shared_environment_variance: float
    residual_variance: float
    h2: float
    c2: float
    beta: FloatArray
    converged: bool
    iterations: int
    objective: Optional[float] = None
    message: str = ""
    se_h2: Optional[float] = None
    se_c2: Optional[float] = None
    se_h2_plus_c2: Optional[float] = None
    covariance_proportions: Optional[FloatArray] = None
    hessian_logratio: Optional[FloatArray] = None
    hessian_condition_number: Optional[float] = None
    hessian_min_eigenvalue: Optional[float] = None
    se_message: str = "not requested"


class ExactProfileML:
    """Exact profile ML using block inverses and a rank x rank Cholesky."""

    def __init__(self, covariance: LowRankACE):
        self.cov = covariance

    def _solve_logdet(self, h2: float, c2: float, rhs: ArrayLike):
        e2 = 1.0 - h2 - c2
        if h2 < 0.0 or c2 < 0.0 or e2 <= 1e-9:
            raise np.linalg.LinAlgError("Variance proportions are outside the simplex.")
        base_diag = e2 + h2 * self.cov.delta
        rinv_rhs = self.cov.families.base_inverse(rhs, base_diag, c2)
        if h2 <= 1e-14 or self.cov.rank == 0:
            return rinv_rhs, self.cov.families.base_logdet(base_diag, c2)

        # R^-1 B is materialized only for the exact backend. Its size is N x r.
        rinv_b = self.cov.families.base_inverse(self.cov.B, base_diag, c2)
        gram = np.eye(self.cov.rank) + h2 * self.cov.bt_matmat(rinv_b)
        cf = cho_factor(gram, lower=True, check_finite=False)
        middle = cho_solve(cf, self.cov.bt_matmat(rinv_rhs), check_finite=False)
        solved = rinv_rhs - h2 * (rinv_b @ middle)
        logdet = self.cov.families.base_logdet(base_diag, c2)
        logdet += 2.0 * np.log(np.diag(cf[0])).sum()
        return solved, float(logdet)

    def fit(
        self,
        y: ArrayLike,
        x: Optional[ArrayLike] = None,
        fix_c2_zero: bool = False,
        compute_standard_errors: bool = False,
        hessian_step: float = 2e-3,
    ) -> FitResult:
        y = np.asarray(y, dtype=np.float64).reshape(-1)
        if y.size != self.cov.n:
            raise ValueError("y must have length N.")
        x = np.ones((self.cov.n, 1)) if x is None else np.asarray(x, dtype=np.float64)
        if x.ndim == 1:
            x = x[:, None]
        rhs = np.column_stack((y, x))

        def unpack(z):
            if fix_c2_zero:
                h2 = 1.0 / (1.0 + np.exp(-float(z[0])))
                return np.array([h2, 0.0, 1.0 - h2])
            ez = np.exp(np.r_[z, 0.0] - np.max(np.r_[z, 0.0]))
            return ez / ez.sum()

        def objective(z):
            h2, c2, _ = unpack(z)
            try:
                solved, logdet = self._solve_logdet(h2, c2, rhs)
                vy, vx = solved[:, 0], solved[:, 1:]
                xtvx = x.T @ vx
                beta = np.linalg.solve(xtvx, x.T @ vy)
                resid = y - x @ beta
                q = float(resid @ (vy - vx @ beta))
                if not np.isfinite(q) or q <= 0:
                    return 1e100
                return 0.5 * (logdet + self.cov.n * np.log(q / self.cov.n))
            except (np.linalg.LinAlgError, ValueError):
                return 1e100

        initial = np.array([np.log(0.4 / 0.6)]) if fix_c2_zero else np.log([0.4, 0.2])
        opt = minimize(objective, initial, method="L-BFGS-B")
        h2, c2, e2 = unpack(opt.x)
        solved, _ = self._solve_logdet(h2, c2, rhs)
        vy, vx = solved[:, 0], solved[:, 1:]
        beta = np.linalg.solve(x.T @ vx, x.T @ vy)
        resid = y - x @ beta
        scale = float(resid @ (vy - vx @ beta)) / self.cov.n
        result = FitResult(
            scale * h2, scale * c2, scale * e2, h2, c2, beta,
            bool(opt.success), int(opt.nit), float(opt.fun), str(opt.message)
        )
        if compute_standard_errors:
            self._add_profile_standard_errors(
                result, opt.x, objective, unpack, fix_c2_zero, hessian_step
            )
        return result

    @staticmethod
    def _central_hessian(function, point: ArrayLike, step: float) -> FloatArray:
        """Central finite-difference Hessian of a scalar function."""
        z = np.asarray(point, dtype=np.float64)
        if step <= 0 or not np.isfinite(step):
            raise ValueError("hessian_step must be a positive finite number.")
        p = z.size
        hessian = np.empty((p, p), dtype=np.float64)
        f0 = float(function(z))
        if not np.isfinite(f0):
            raise np.linalg.LinAlgError("Objective is not finite at the optimum.")
        eye = np.eye(p) * step
        for i in range(p):
            fp = float(function(z + eye[i]))
            fm = float(function(z - eye[i]))
            hessian[i, i] = (fp - 2.0 * f0 + fm) / (step * step)
            for j in range(i):
                fpp = float(function(z + eye[i] + eye[j]))
                fpm = float(function(z + eye[i] - eye[j]))
                fmp = float(function(z - eye[i] + eye[j]))
                fmm = float(function(z - eye[i] - eye[j]))
                value = (fpp - fpm - fmp + fmm) / (4.0 * step * step)
                hessian[i, j] = value
                hessian[j, i] = value
        return 0.5 * (hessian + hessian.T)

    @staticmethod
    def _proportion_jacobian(h2: float, c2: float, fix_c2_zero: bool) -> FloatArray:
        """Jacobian of (h2, c2, h2+c2) with respect to log ratios."""
        if fix_c2_zero:
            derivative = h2 * (1.0 - h2)
            return np.array([[derivative], [0.0], [derivative]], dtype=np.float64)
        return np.array([
            [h2 * (1.0 - h2), -h2 * c2],
            [-h2 * c2, c2 * (1.0 - c2)],
            [h2 * (1.0 - h2 - c2), c2 * (1.0 - h2 - c2)],
        ], dtype=np.float64)

    def _add_profile_standard_errors(
        self, result: FitResult, optimum: ArrayLike, objective, unpack,
        fix_c2_zero: bool, hessian_step: float,
    ) -> None:
        """Add observed-information/delta-method SEs to a fitted result."""
        try:
            hessian = self._central_hessian(objective, optimum, hessian_step)
            eigenvalues = np.linalg.eigvalsh(hessian)
            min_eigenvalue = float(eigenvalues.min())
            condition = float(np.linalg.cond(hessian))
            if min_eigenvalue <= 0 or not np.isfinite(condition):
                raise np.linalg.LinAlgError(
                    "Profile Hessian is not positive definite; Wald SEs are unreliable."
                )
            covariance_logratio = np.linalg.inv(hessian)
            h2, c2, _ = unpack(np.asarray(optimum, dtype=np.float64))
            jacobian = self._proportion_jacobian(float(h2), float(c2), fix_c2_zero)
            covariance_proportions = jacobian @ covariance_logratio @ jacobian.T
            variances = np.diag(covariance_proportions)
            if np.any(variances < -1e-10):
                raise np.linalg.LinAlgError("Delta-method covariance has a negative diagonal.")
            standard_errors = np.sqrt(np.maximum(variances, 0.0))
            result.se_h2 = float(standard_errors[0])
            result.se_c2 = None if fix_c2_zero else float(standard_errors[1])
            result.se_h2_plus_c2 = float(standard_errors[2])
            result.covariance_proportions = covariance_proportions
            result.hessian_logratio = hessian
            result.hessian_condition_number = condition
            result.hessian_min_eigenvalue = min_eigenvalue
            result.se_message = "observed profile-likelihood Hessian; delta method"
        except (ValueError, np.linalg.LinAlgError, FloatingPointError) as error:
            result.se_message = f"standard errors unavailable: {error}"


class MatrixFreeAIREML:
    """Matrix-free stochastic AI-REML for large N.

    The trace terms use fixed Rademacher probes. Increase n_trace_probes for
    final inference; 20-50 is a useful range, while 8-12 is suitable for a
    preliminary run. This estimator never factors an r x r matrix.
    """

    def __init__(
        self,
        covariance: LowRankACE,
        n_trace_probes: int = 20,
        cg_tolerance: float = 1e-5,
        cg_maxiter: int = 500,
        random_state: int = 1,
    ):
        self.cov = covariance
        self.n_trace_probes = int(n_trace_probes)
        self.cg_tolerance = float(cg_tolerance)
        self.cg_maxiter = int(cg_maxiter)
        self.rng = np.random.default_rng(random_state)

    def _solver(self, theta: FloatArray):
        n = self.cov.n
        vop = LinearOperator((n, n), matvec=lambda z: self.cov.covariance(theta, z))
        base_diag = theta[2] + theta[0] * self.cov.delta
        mop = LinearOperator(
            (n, n),
            matvec=lambda z: self.cov.families.base_inverse(z, base_diag, theta[1]),
        )

        def solve(rhs):
            a, was_vector = _as_2d(rhs, n)
            out = np.empty_like(a)
            for j in range(a.shape[1]):
                out[:, j], info = cg(
                    vop, a[:, j], M=mop, rtol=self.cg_tolerance,
                    atol=0.0, maxiter=self.cg_maxiter
                )
                if info != 0:
                    raise RuntimeError(f"PCG did not converge (info={info}).")
            return _restore_shape(out, was_vector)

        return solve

    def fit(
        self,
        y: ArrayLike,
        x: Optional[ArrayLike] = None,
        initial: Optional[ArrayLike] = None,
        max_iterations: int = 30,
        tolerance: float = 2e-3,
    ) -> FitResult:
        y = np.asarray(y, dtype=np.float64).reshape(-1)
        if y.size != self.cov.n:
            raise ValueError("y must have length N.")
        x = np.ones((self.cov.n, 1)) if x is None else np.asarray(x, dtype=np.float64)
        if x.ndim == 1:
            x = x[:, None]
        if x.shape[0] != self.cov.n:
            raise ValueError("X must have N rows.")
        total = float(np.var(y, ddof=max(1, x.shape[1])))
        theta = np.asarray(
            initial if initial is not None else [0.4 * total, 0.1 * total, 0.5 * total],
            dtype=np.float64,
        )
        if theta.shape != (3,) or np.any(theta <= 0):
            raise ValueError("initial must contain three positive variances.")
        probes = self.rng.choice(
            np.array([-1.0, 1.0]), size=(self.cov.n, self.n_trace_probes)
        )
        converged = False
        message = "Maximum iterations reached."
        beta = np.zeros(x.shape[1])

        for iteration in range(1, max_iterations + 1):
            solve = self._solver(theta)
            vinvx = solve(x)
            xtvinvx_cf = cho_factor(x.T @ vinvx, lower=True, check_finite=False)

            def apply_p(rhs):
                vinvr = solve(rhs)
                a, was_vector = _as_2d(vinvr, self.cov.n)
                adjustment = vinvx @ cho_solve(
                    xtvinvx_cf, x.T @ a, check_finite=False
                )
                return _restore_shape(a - adjustment, was_vector)

            q = apply_p(y)
            beta = cho_solve(xtvinvx_cf, x.T @ solve(y), check_finite=False)
            pz = apply_p(probes)
            traces = np.empty(3)
            quadratic = np.empty(3)
            w = []
            for k in range(3):
                kz = self.cov.component(k, probes)
                traces[k] = np.einsum("ij,ij->", pz, kz) / self.n_trace_probes
                wk = self.cov.component(k, q)
                w.append(wk)
                quadratic[k] = q @ wk
            score = 0.5 * (quadratic - traces)
            pw = [apply_p(wk) for wk in w]
            ai = np.empty((3, 3))
            for i in range(3):
                for j in range(i, 3):
                    ai[i, j] = ai[j, i] = 0.5 * (w[i] @ pw[j])
            ridge = max(np.trace(ai), 1.0) * 1e-10
            step = np.linalg.solve(ai + ridge * np.eye(3), score)
            step_scale = 1.0
            while np.any(theta + step_scale * step <= theta * 1e-5):
                step_scale *= 0.5
                if step_scale < 2.0 ** -20:
                    break
            new_theta = np.maximum(theta + step_scale * step, theta * 1e-5)
            relative_change = np.max(np.abs(new_theta - theta) / np.maximum(theta, 1e-12))
            theta = new_theta
            if relative_change < tolerance:
                converged = True
                message = "Relative variance-component change met tolerance."
                break

        total_var = float(theta.sum())
        return FitResult(
            theta[0], theta[1], theta[2], theta[0] / total_var,
            theta[1] / total_var, beta, converged, iteration, None, message
        )


def nystrom_from_cross_grm(
    cross_grm: ArrayLike,
    knot_grm: ArrayLike,
    grm_diagonal: ArrayLike,
    jitter: float = 1e-8,
) -> tuple[FloatArray, FloatArray]:
    """Convert A[:, knots], A[knots, knots], diag(A) to K=B B'+diag(delta)."""
    c = np.asarray(cross_grm, dtype=np.float64)
    w = np.asarray(knot_grm, dtype=np.float64)
    diag = np.asarray(grm_diagonal, dtype=np.float64)
    if c.ndim != 2 or w.shape != (c.shape[1], c.shape[1]) or diag.shape != (c.shape[0],):
        raise ValueError("Incompatible cross-GRM, knot-GRM, or diagonal dimensions.")
    ridge = jitter * max(float(np.trace(w) / max(w.shape[0], 1)), 1.0)
    chol = np.linalg.cholesky(w + ridge * np.eye(w.shape[0]))
    # B = C L^{-T}, evaluated by solving L B' = C'.
    factor = np.linalg.solve(chol, c.T).T
    explained_diag = np.einsum("ij,ij->i", factor, factor)
    delta = np.maximum(diag - explained_diag, 0.0)
    return factor, delta
