cmake_minimum_required(VERSION 3.10)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS YES CACHE BOOL "Export all symbols")

add_subdirectory(extern)

project(svd VERSION 1.0.0 DESCRIPTION "Implementations of SVD")

SET(SOURCE_FILES
  src/GolubKahanSVD.cpp
  #src/iterative_refinement.cpp
  src/jacobi.cpp
  src/generate_svd.cpp
  #src/testing.cpp
)
SET(HEADER_FILES
  src/generate_svd.h
  src/jacobi.h
)

add_library(svd ${SOURCE_FILES} ${HEADER_FILES}) 
target_compile_features(svd PRIVATE cxx_std_20)
target_link_libraries(svd eigen)

install(
  TARGETS svd
  LIBRARY DESTINATION lib
)
install(
  FILES ${HEADER_FILES} 
  DESTINATION include/library_name
)
