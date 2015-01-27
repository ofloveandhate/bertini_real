




#we're building bertini_real, version 0.1.1, and the corresponding email is dan brake's
AC_INIT([bertini_real], [0.1.1], [danielthebrake@gmail.com],[bertini_real], [http://www.danielthebrake.org/])


# Force autoconf to be at least this version number:
AC_PREREQ([2.68])


#
AC_CONFIG_AUX_DIR([config])

# turn on the keeping of produced objects in their folders.  this is for non-recursive make
# and autotools
# see Recursive Make Considered Harmful, and any number of demos.
AM_INIT_AUTOMAKE([1.13 subdir-objects])

#the language for BertiniReal is C++
AC_LANG([C++])

#use the local m4 folder.
ACLOCAL_AMFLAGS = -I m4

#another directive to use the m4 folder
AC_CONFIG_MACRO_DIR([m4])

#the only produced file will be a single Makefile.
AC_CONFIG_FILES([Makefile])

#find the CC compiler
AC_PROG_CXX

#set up for building libraries
AC_PROG_LIBTOOL

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

#enable the creation of shared libraries
AC_ENABLE_SHARED

#enable the creation of static libraries
AC_ENABLE_STATIC



#this calls a file in the m4/ directory, which sets up the MPI wrapper stuffs
LX_FIND_MPI

AS_IF([test "$have_CXX_mpi" = "yes"],
	AC_MSG_NOTICE([mpi c++ appears to work correctly]),
	AC_MSG_ERROR([it appears you do not have mpi installed])
	)



# the form of the following commands --
# AC_SEARCH_LIBS(function, libraries-list, action-if-found, action-if-not-found, extra-libraries)

# find mpfr
AC_SEARCH_LIBS([mpfr_get_version],[mpfr], [],[
	AC_MSG_ERROR([unable to find mpfr])
])

#find gmp
AC_SEARCH_LIBS([__gmpz_init],[gmp], [],[
	AC_MSG_ERROR([unable to find gmp])
])

#find bertini
AC_SEARCH_LIBS([head_zero_dim_track_d],[bertini-parallel], [],[
	AC_MSG_ERROR([unable to find libbertini-parallel._])
])

#find the math library
AC_SEARCH_LIBS([cos], [m], [], [
  AC_MSG_ERROR([unable to find the cos() function])
  ])



AX_BOOST_BASE([1.53],, [AC_MSG_ERROR([bertini_real needs Boost, but it was not found in your system])])
AX_BOOST_SYSTEM
AX_BOOST_FILESYSTEM
AX_BOOST_CHRONO
AX_BOOST_REGEX
AX_BOOST_TIMER





#wrap it up.
AC_OUTPUT







