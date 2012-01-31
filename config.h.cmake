#ifndef CONFIG_H
#define CONFIG_H

#define DUNE_MINIMAL_DEBUG_LEVEL 4
#cmakedefine HAVE_BOOST 1
#cmakedefine HAVE_DUNE 1
#cmakedefine HAVE_DUNE_GRID 1
#cmakedefine HAVE_DUNE_DISC 1
#cmakedefine HAVE_DUNE_FEM 1
#cmakedefine HAVE_DUNE_ISTL 1
#cmakedefine HAVE_DUNE_LOCALFUNCTIONS 1
#cmakedefine HAVE_DUNE_PDELAB 1

#ifdef ENABLE_MPI
#cmakedefine HAVE_MPI 1
#endif

#ifdef ENABLE_UG
#cmakedefine HAVE_UG 1
#ifdef ENABLE_MPI
#define ModelP
#endif
#endif

#ifdef ENABLE_ALUGRID
#cmakedefine HAVE_ALUGRID 1
#endif

#ifdef ENABLE_METIS
#cmakedefine HAVE_METIS 1
#endif

#ifdef ENABLE_ALBERTA
#cmakedefine HAVE_ALBERTA 1
#endif

#cmakedefine PROJECT_NAME             "${PROJECT_NAME}"
#cmakedefine PROJECT_VERSION          "${PROJECT_VERSION}"
#cmakedefine PROJECT_MAINTAINER       "${PROJECT_MAINTAINER}"
#cmakedefine PROJECT_MAINTAINER_EMAIL "${PROJECT_MAINTAINER_EMAIL}"

#cmakedefine HAVE_SUPERLU ENABLE_SUPERLU

/* tr1/array. */
//#cmakedefine HAVE_TR1_ARRAY 1
#cmakedefine HAVE_NULLPTR 1
#cmakedefine HAVE_STATIC_ASSERT 1
#cmakedefine HAVE_VARIADIC_TEMPLATES 1
#cmakedefine HAVE_VARIADIC_CONSTRUCTOR_SFINAE 1
#cmakedefine HAVE_RVALUE_REFERENCES 1
#cmakedefine HAVE_MALLOC_H 1
#cmakedefine HAVE_VALGRIND 1

#define DUMUX_DEPRECATED __attribute__((deprecated))
#define DUMUX_DEPRECATED_MSG(msg) __attribute__((deprecated(msg)))

#define DUNE_DEPRECATED __attribute__((deprecated))
#define DUNE_DEPRECATED_MSG(msg) __attribute__((deprecated(msg)))

#endif // CONFIG_H
