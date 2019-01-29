file(GLOB_RECURSE GAMMA_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_executable(gamma ${GAMMA_SOURCES})

message(${PUGIXML_INCLUDE_DIR})
target_include_directories(gamma PUBLIC ${EIGEN3_INCLUDE_DIR})

target_link_libraries(gamma ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(gamma ${MPI_LIBRARIES})
   target_link_libraries(gamma ${BLAS_LIBRARIES})
target_link_libraries(gamma Faddeeva)
target_link_libraries(gamma libpugixml)
target_link_libraries(gamma libctf)
target_link_libraries(gamma libint)
target_link_libraries(gamma libecpint)

if(MPI_COMPILE_FLAGS)
  set_target_properties(gamma PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(gamma PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

install(TARGETS gamma DESTINATION bin)
