# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

# tests that should build and run successfully
link_libraries(Dune::Functions)

dune_add_test(SOURCES containerfactorytest.cc LABELS quick)

dune_add_test(SOURCES istlvectorbackendtest.cc LABELS quick)

dune_add_test(SOURCES istlvectorfactorytest.cc LABELS quick)
