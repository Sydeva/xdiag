#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 The xdiag Contributors
"""Smoke test for genlat_triangular_tonly.py."""

import subprocess
import tempfile
from pathlib import Path

import tomllib


def main():
    root = Path(__file__).resolve().parent
    script = root / "genlat_triangular_tonly.py"

    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "triangular.test.toml"
        cmd = [
            "python3",
            str(script),
            "4",
            "2",
            "--verify",
            "--toml",
            str(out),
        ]
        subprocess.run(cmd, check=True, cwd=root)

        assert out.exists() and out.stat().st_size > 0
        data = tomllib.loads(out.read_text(encoding="utf-8"))

        assert "Translation_a1" in data
        assert "Translation_2a2" in data
        assert "Translation_a2" not in data

        syms = data["Symmetries"]
        irreps = [k for k in data.keys() if k.startswith("k") and k[1:].isdigit()]
        assert len(irreps) == len(syms)

        for name in irreps:
            section = data[name]
            assert "characters" in section
            assert "momentum" in section
            assert "allowed_symmetries" not in section
            assert len(section["characters"]) == len(syms)

        nnn_pairs = set()
        ssc_terms = []
        for term in data["Interactions"]:
            if term[0] != "J2":
                if len(term) == 5 and term[1] == "ScalarChirality":
                    ssc_terms.append(term)
                continue
            i, j = int(term[2]), int(term[3])
            pair = (i, j) if i < j else (j, i)
            assert pair not in nnn_pairs
            nnn_pairs.add(pair)

        # ScalarChirality terms should be emitted in oriented triplets.
        assert len(ssc_terms) > 0
        assert len(ssc_terms) % 3 == 0
        for term in ssc_terms:
            assert term[1] == "ScalarChirality"
            i, j, k = int(term[2]), int(term[3]), int(term[4])
            assert i != j and j != k and i != k

    print("Smoke test passed.")


if __name__ == "__main__":
    main()

