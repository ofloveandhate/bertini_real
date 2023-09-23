




#we're building bertini_real, version 1.6.1, and the corresponding email is silviana's
AC_INIT([bertini_real],[1.6.1],[silviana.amethyst@gmail.com],[bertini_real],[https://bertinireal.com/])


# Force autoconf to be at least this version number:
AC_PREREQ([2.71])


#
AC_CONFIG_AUX_DIR([config])

# turn on the keeping of produced objects in their folders.  this is for non-recursive make
# and autotools
# see Recursive Make Considered Harmful, and any number of demos.
AM_INIT_AUTOMAKE([1.13 subdir-objects])

#the language for BertiniReal is C++
AC_LANG([C++])

#another directive to use the m4 folder
AC_CONFIG_MACRO_DIR([m4])

#the only produced file will be a single Makefile.
AC_CONFIG_FILES([Makefile])

#find the CC compiler
AC_PROG_CXX

#this calls a file in the m4/ directory, which sets up the MPI wrapper stuffs

# If --with-mpi=auto is used, try to find MPI, but use standard C compiler if it is not found.
# If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
# If --with-mpi=no is used, use a standard C compiler instead.
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
    [compile with MPI (parallelization) support. If none is found,
    MPI is not used. Default: auto])
],,[with_mpi=auto])

AX_PROG_CXX_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
  use_mpi=no
  if test x"$with_mpi" = xyes; then
    AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
  else
    AC_MSG_FAILURE([No MPI compiler found, won't use MPI.])
  fi
])




#enable the creation of shared libraries
AC_ENABLE_SHARED

#enable the creation of static libraries
AC_ENABLE_STATIC

#set up for building libraries
LT_INIT

#find flex / lex
AC_PROG_LEX

#find bison / yacc
AC_PROG_YACC

#find the linker
AC_PROG_LN_S

#find the command for making directories with their parents
AC_PROG_MKDIR_P

#fire up libtool
LT_INIT


AX_CXX_COMPILE_STDCXX_11






# the form of the following commands --
# AC_SEARCH_LIBS(function, libraries-list, action-if-found, action-if-not-found, extra-libraries)

#find gmp
AC_SEARCH_LIBS([__gmpz_init],[gmp], [],[
	AC_MSG_ERROR([unable to find gmp])
])

# find mpfr
AC_SEARCH_LIBS([mpfr_get_version],[mpfr], [],[
	AC_MSG_ERROR([unable to find mpfr])
])



#find bertini
AC_SEARCH_LIBS([head_zero_dim_track_d],[bertini-parallel], [],
	[AC_MSG_ERROR([unable to find libbertini-parallel.  please build and install bertini1 from source.  bertini.nd.edu])],
	$MPI_CXXLDFLAGS)

#find the math library
AC_SEARCH_LIBS([cos], [m], [], [
  AC_MSG_ERROR([unable to find the cos() function])
  ])



AX_BOOST_BASE([1.50],, [AC_MSG_ERROR([bertini_real needs Boost, but it was not found in your system])])

AX_BOOST_FILESYSTEM

AX_BOOST_REGEX
AX_BOOST_CHRONO
AX_BOOST_TIMER

AX_BOOST_SYSTEM


AM_CONFIG_HEADER(config.h)


#wrap it up.
AC_OUTPUT
