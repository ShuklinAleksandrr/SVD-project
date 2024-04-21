cmake_minimum_required(VERSION 3.10)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS YES CACHE BOOL "Export all symbols")

add_subdirectory(extern)

project(svd VERSION 1.0.0 DESCRIPTION "Implementations of SVD")

SET(SOURCE_FILES
  src/GolubKahanSVD.cpp
  src/iterative_refinement.cpp
  src/jacobi.cpp
  src/generate_svd.cpp
)
SET(HEADER_FILES
  src/iterative_refinement.h
  src/generate_svd.h
  src/jacobi.h
  src/dqds.h
)

add_library(svd ${SOURCE_FILES} ${HEADER_FILES}) 
target_compile_features(svd PRIVATE cxx_std_20)
target_link_libraries(svd eigen)
target_link_libraries(svd lapacke)

install(
  TARGETS svd
  LIBRARY DESTINATION lib
)
install(
  FILES ${HEADER_FILES} 
  DESTINATION include/library_name
)

add_executable (svd_test src/testing.cpp)
target_compile_features(svd_test PRIVATE cxx_std_20)
target_link_libraries(svd_test eigen)
target_link_libraries(svd_test lapacke)
