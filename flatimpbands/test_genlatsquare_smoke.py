#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 The xdiag Contributors

import tempfile
from pathlib import Path

from genlatsquare import (
    compute_irrep_sections,
    compute_translation_perm,
    generate_all_translations,
    verify_data,
    write_toml,
)


def main():
    lx, ly = 3, 2

    ta1 = compute_translation_perm(lx, ly, 1, 0)
    ta2 = compute_translation_perm(lx, ly, 0, 1)
    assert len(ta1) == lx * ly
    assert len(ta2) == lx * ly

    symmetries, displacements = generate_all_translations(lx, ly)
    sections = compute_irrep_sections(lx, ly, displacements)
    verify_data(lx, ly, symmetries, sections)
    assert len(sections) == lx * ly
    assert sections[0]["name"] == "k1"

    with tempfile.TemporaryDirectory() as tmpdir:
        out = Path(tmpdir) / "square_test.toml"
        write_toml(str(out), lx, ly)
        text = out.read_text(encoding="utf-8")
        assert "Coordinates = [" in text
        assert "Symmetries = [" in text
        assert "[k1]" in text
        assert "allowed_symmetries" not in text

    print("test_genlatsquare_smoke: OK")


if __name__ == "__main__":
    main()

