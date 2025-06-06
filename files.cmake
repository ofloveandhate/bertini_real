set(common_headers
	include/bertini1/bertini_extensions.hpp
	include/bertini1/bertini_headers.hpp

	include/cells/cell.hpp
	include/cells/edge.hpp
	include/cells/face.hpp
	include/cells/vertex.hpp

	include/containers/holders.hpp
	include/containers/vertex_set.hpp

	include/decompositions/checkSelfConjugate.hpp
	include/decompositions/curve.hpp
	include/decompositions/decomposition.hpp
	include/decompositions/surface.hpp

	include/io/color.hpp
	include/io/fileops.hpp
	include/io/partitionParse.h

	include/nag/solvers/midpoint.hpp
	include/nag/solvers/multilintolin.hpp
	include/nag/solvers/nullspace.hpp
	include/nag/solvers/postProcessing.hpp
	include/nag/solvers/solver.hpp
	include/nag/solvers/sphere_intersection.hpp

	include/nag/nid.hpp
	include/nag/system_randomizer.hpp
	include/nag/witness_set.hpp

	include/symbolics/derivative_systems.hpp
	include/symbolics/isosingular.hpp
	include/symbolics/nullspace.hpp
	include/symbolics/slicing.hpp
	include/symbolics/sphere_intersection.hpp

	include/double_odometer.hpp
	include/forward_declarations.hpp
	include/limbo.hpp
	include/parallelism.hpp
	include/programConfiguration.hpp
)


set(common_src
	src/bertini1/bertini_extensions.cpp

	src/cells/face.cpp
	src/cells/vertex.cpp

	src/containers/holders.cpp
	src/containers/vertex_set.cpp

	src/decompositions/checkSelfConjugate.cpp
	src/decompositions/curve.cpp
	src/decompositions/decomposition.cpp
	src/decompositions/surface.cpp

	src/io/color.cpp
	src/io/fileops.cpp

	src/nag/solvers/midpoint.cpp
	src/nag/solvers/multilintolin.cpp
	src/nag/solvers/nullspace.cpp
	src/nag/solvers/postProcessing.cpp
	src/nag/solvers/solver.cpp
	src/nag/solvers/sphere_intersection.cpp

	src/nag/nid.cpp
	src/nag/system_randomizer.cpp
	src/nag/witness_set.cpp

	src/symbolics/derivative_systems.cpp
	src/symbolics/isosingular.cpp
	src/symbolics/nullspace.cpp
	src/symbolics/slicing.cpp
	src/symbolics/sphere_intersection.cpp

	src/parallelism.cpp
	src/programConfiguration.cpp
)



set(sampler_headers
	include/sampler.hpp
)

set(br_headers
	include/bertini_real.hpp
)



set(sampler_src
	src/sampler/curve_methods.cpp
	src/sampler/surface_methods.cpp
	src/sampler/sampler.cpp
)

set(br_src
	src/bertini_real.cpp
)