CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

# define executable
add_executable(two_d src/main.cc
		    )


target_link_libraries(two_d LF::lf.base
			    LF::lf.mesh.test_utils
			    LF::lf.mesh.hybrid2d
			    Boost::filesystem
			    Boost::system
		      	    Boost::program_options
			    LF::lf.uscalfe
			    LF::lf.quad
			    LF::lf.geometry
			    LF::lf.assemble
			    LF::lf.io)

