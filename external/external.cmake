add_custom_target(External)

add_subdirectory(external/Faddeeva)
include_directories(${CMAKE_SOURCE_DIR}/external/Faddeeva)
set(FADDEEVA_LIBRARY ${CMAKE_BINARY_DIR}/external/Faddeeva/libFaddeeva${CMAKE_STATIC_LIBRARY_SUFFIX})

set(EXTERNAL_SOURCE_DIR ${PROJECT_SOURCE_DIR}/external)
set(EXTERNAL_BUILD_DIR ${PROJECT_BINARY_DIR}/external/build)

#
# add pugixml
#
find_library(PUGIXML_LIBRARY pugixml
			HINTS "${LIBPATH}" "${LIBPATH}/lib" "/usr/lib" "/usr/local/lib")
find_path(PUGIXML_INCLUDE_DIR pugixml.hpp
			HINTS "${INCPATH}/include" "/usr/include" "/usr/local/include")

if((NOT PUGIXML_LIBRARY) OR (NOT PUGIXML_INCLUDE_DIR))
	message("Unable to find pugixml, cloning...")
	
	ExternalProject_Add(pugixml_external
			PREFIX ${EXTERNAL_BUILD_DIR}/pugixml
			GIT_REPOSITORY https://github.com/zeux/pugixml
			CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
			UPDATE_COMMAND ""
		)
	
	set(PUGIXML_LIBRARY ${EXTERNAL_BUILD_DIR}/pugixml/src/pugixml_external-build/libpugixml${CMAKE_STATIC_LIBRARY_SUFFIX})
	set(PUGIXML_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/build/pugixml/src/pugixml_external/src)
	
	add_dependencies(External pugixml_external)
	add_library(libpugixml STATIC IMPORTED)
	set_target_properties(libpugixml PROPERTIES IMPORTED_LOCATION ${PUGIXML_LIBRARY})
	
else()
	message("Found pugixml")
	message("${PUGIXML_LIBRARY}")
	
	add_library(libpugixml INTERFACE)
	target_include_directories(libpugixml INTERFACE ${PUGIXML_INCLUDE_DIR})
	target_link_libraries(libpugixml INTERFACE ${PUGIXML_LIBRARY})
endif()

include_directories(${PUGIXML_INCLUDE_DIR})

#
# always pull and build ctf as it has no sensible install setup
#

set(CTF_BLAS_LIBS "-lblas")
set(CTF_INCLUDES "")
set(CTF_DEFS "")

	ExternalProject_Add(ctf_external
			PREFIX ${EXTERNAL_BUILD_DIR}/ctf
			GIT_REPOSITORY https://github.com/solomonik/ctf
			GIT_TAG 2b5d0b8
			BUILD_IN_SOURCE 1
			UPDATE_COMMAND ""
			CONFIGURE_COMMAND /bin/bash ./configure
			"CXXFLAGS=${CMAKE_CXX_FLAGS} -DOMP_OFF"
			"WARNFLAGS=-Wall -Wno-format"
			"--blas=${CTF_BLAS_LIBS}"
			"INCLUDES=${CTF_INCLUDES}"			
			BUILD_COMMAND $(MAKE) SHELL=/bin/bash
			INSTALL_COMMAND ""
		)
add_dependencies(External ctf_external)	
set(CTF_LIBRARY ${EXTERNAL_BUILD_DIR}/ctf/src/ctf_external/lib/libctf${CMAKE_STATIC_LIBRARY_SUFFIX})
set(CTF_INCLUDE_DIR ${EXTERNAL_BUILD_DIR}/ctf/src/ctf_external/include)

add_library(libctf STATIC IMPORTED)
set_target_properties(libctf PROPERTIES IMPORTED_LOCATION ${CTF_LIBRARY})
include_directories(${CTF_INCLUDE_DIR})


#
# find Libint
#
find_library(LIBINT2_LIBRARY int2
			HINTS ${LIBPATH}/lib ${LIBPATH}/libint/2.4.0-beta.2/lib "/usr/local/libint/2.4.0-beta.2/lib")
find_path(LIBINT2_INCLUDE_DIR libint2.h
			HINTS ${INCPATH}/include ${INCPATH}/libint/2.4.0-beta.2/include "/usr/local/libint/2.4.0-beta.2/include")

if ((NOT LIBINT2_LIBRARY) OR (NOT LIBINT2_INCLUDE_DIR))	
	message("Unable to find Libint2, untarring...")
	
	ExternalProject_Add(libint_external
			PREFIX ${EXTERNAL_BUILD_DIR}/libint
			URL ${CMAKE_SOURCE_DIR}/external/libint/libint.tgz
			BUILD_IN_SOURCE 1
			CONFIGURE_COMMAND /bin/bash ./configure 
			--prefix=${EXTERNAL_BUILD_DIR}/libint
			"CFLAGS=${CMAKE_C_FLAGS}"
			"CXXFLAGS=${CMAKE_CXX_FLAGS}"
			"CPPFLAGS=-I${EIGEN3_INCLUDE_DIR}"
			BUILD_COMMAND $(MAKE) SHELL=/bin/bash
			INSTALL_COMMAND make install SHELL=/bin/bash
		)
		
	add_dependencies(External libint_external)
	set(LIBINT2_LIBRARY ${EXTERNAL_BUILD_DIR}/libint/src/libint_external/lib/libint2${CMAKE_STATIC_LIBRARY_SUFFIX})
	set(LIBINT2_INCLUDE_DIR ${EXTERNAL_BUILD_DIR}/libint/src/libint_external/include)
else()
	message("Found Libint2")
	message("${LIBINT2_LIBRARY}")	
endif()

add_library(libint STATIC IMPORTED)
set_target_properties(libint PROPERTIES IMPORTED_LOCATION ${LIBINT2_LIBRARY})
include_directories(${LIBINT2_INCLUDE_DIR})
include_directories(${LIBINT2_INCLUDE_DIR}/libint2)

#
# find Libecpint
#
find_library(LIBECPINT_LIBRARY ecpint
			HINTS "${LIBPATH}/lib" "/usr/lib" "/usr/local/lib")
find_path(LIBECPINT_INCLUDE_DIR libecpint.hpp
			HINTS "/usr/local/include" "/usr/include" "${INCPATH}/include")

if ((NOT LIBECPINT_LIBRARY) OR (NOT LIBECPINT_INCLUDE_DIR))
	message("Unable to find Libecpint, cloning…")

	ExternalProject_Add(libecpint_external
			    PREFIX ${EXTERNAL_BUILD_DIR}
			    GIT_REPOSITORY https://github.com/robashaw/libecpint
			    CMAKE_ARGS "-DLIBECPINT_MAX_L=6 -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
			    UPDATE_COMMAND ""
			)
	set(LIBECPINT_LIBRARY ${EXTERNAL_BUILD_DIR}/libecpint/src/libecpint_external-build/libecpint${CMAKE_STATIC_LIBRARY_SUFFIX})
	set(LIBECPINT_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/libecpint/src/libecpint_external/include/libecpint)

	add_dependencies(External libecpint_external)	
	add_library(libecpint STATIC IMPORTED)
	set_target_properties(libecpint PROPERTIES IMPORTED_LOCATION ${LIBECPINT_LIBRARY})
	
else()
	message("Found Libecpint")
	message("${LIBECPINT_LIBRARY}")
	

	add_library(libecpint INTERFACE)
	target_include_directories(libecpint INTERFACE ${LIBECPINT_INCLUDE_DIR})
	target_link_libraries(libecpint INTERFACE ${LIBECPINT_LIBRARY})
endif()