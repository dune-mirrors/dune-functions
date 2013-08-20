// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FINAL_HH
#define DUNE_FUNCTIONS_FINAL_HH

/** \file
 * \brief Definition of the DUNE_FINAL macro, which encapsulates
 * the 'final' keyword from C++11
 */

#if ! HAS_KEYWORD_FINAL
#define DUNE_FINAL
#else
#define DUNE_FINAL final
#endif

#endif
