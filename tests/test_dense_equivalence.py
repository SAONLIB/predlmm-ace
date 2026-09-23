import numpy as np
import tempfile
from pathlib import Path

from predlmm_ace.model import ExactProfileML, LowRankACE, nystrom_from_cross_grm
from predlmm_ace.grm_to_nystrom import build_from_plink_grm


def test_implicit_components_and_block_inverse_match_dense():
    rng = np.random.default_rng(4)
    n, rank = 18, 5
    b = rng.normal(size=(n, rank)) / np.sqrt(rank)
    delta = rng.uniform(0.05, 0.3, n)
    labels = np.repeat(np.arange(6), 3)
    cov = LowRankACE(b, delta, labels, row_chunk_size=7)
    pi = (labels[:, None] == labels[None, :]).astype(float)
    kg = b @ b.T + np.diag(delta)
    x = rng.normal(size=(n, 3))
    np.testing.assert_allclose(cov.genetic(x), kg @ x, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(cov.families.matmat(x), pi @ x, rtol=1e-12, atol=1e-12)

    diagonal, c = rng.uniform(0.2, 1.0, n), 0.17
    base = np.diag(diagonal) + c * pi
    got = cov.families.base_inverse(x, diagonal, c)
    np.testing.assert_allclose(got, np.linalg.solve(base, x), rtol=1e-11, atol=1e-11)
    sign, logdet = np.linalg.slogdet(base)
    assert sign > 0
    np.testing.assert_allclose(cov.families.base_logdet(diagonal, c), logdet, rtol=1e-12)

    h2, c2 = 0.31, 0.14
    v = h2 * kg + c2 * pi + (1.0 - h2 - c2) * np.eye(n)
    exact_solve, exact_logdet = ExactProfileML(cov)._solve_logdet(h2, c2, x)
    np.testing.assert_allclose(exact_solve, np.linalg.solve(v, x), rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(exact_logdet, np.linalg.slogdet(v)[1], rtol=1e-10, atol=1e-10)


def test_nystrom_conversion_reconstructs_diagonal():
    rng = np.random.default_rng(7)
    z = rng.normal(size=(15, 8))
    a = z @ z.T / z.shape[1] + 0.02 * np.eye(15)
    knots = np.array([0, 3, 5, 9])
    b, delta = nystrom_from_cross_grm(a[:, knots], a[np.ix_(knots, knots)], np.diag(a))
    np.testing.assert_allclose(np.einsum("ij,ij->i", b, b) + delta, np.diag(a), atol=2e-7)


def test_plink_triangle_reader_and_ae_fit():
    rng = np.random.default_rng(12)
    n = 24
    z = rng.normal(size=(n, 10))
    a = z @ z.T / 10 + 0.15 * np.eye(n)
    knots = [1, 5, 9, 16, 21]
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        prefix = root / "synthetic"
        with Path(str(prefix) + ".grm.bin").open("wb") as stream:
            for i in range(n):
                np.asarray(a[i, : i + 1], dtype="<f4").tofile(stream)
        with Path(str(prefix) + ".grm.id").open("w") as stream:
            for i in range(n):
                stream.write(f"F{i // 2}\tI{i}\n")
        knot_file = root / "knots.tsv"
        with knot_file.open("w") as stream:
            stream.write("FID\tIID\n")
            for i in knots:
                stream.write(f"F{i // 2}\tI{i}\n")
        out = root / "factor"
        build_from_plink_grm(str(prefix), str(knot_file), str(out), row_chunk_size=7)
        cross = np.load(out / "cross_grm.npy")
        np.testing.assert_allclose(cross, a[:, knots], rtol=2e-7, atol=2e-7)
        factor = np.load(out / "nystrom_factor.npy", mmap_mode="r")
        delta = np.load(out / "diagonal_correction.npy")
        np.testing.assert_allclose(
            np.einsum("ij,ij->i", factor.astype(float), factor.astype(float)) + delta,
            np.diag(a), rtol=1e-7, atol=1e-7
        )
        labels = np.repeat(np.arange(n // 2), 2)
        covariance = LowRankACE(factor, delta, labels)
        y = rng.normal(size=n)
        result = ExactProfileML(covariance).fit(y, np.ones((n, 1)), fix_c2_zero=True)
        assert result.converged
        assert result.c2 == 0.0


if __name__ == "__main__":
    test_implicit_components_and_block_inverse_match_dense()
    test_nystrom_conversion_reconstructs_diagonal()
    test_plink_triangle_reader_and_ae_fit()
    print("All dense-equivalence tests passed.")
