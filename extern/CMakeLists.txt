cmake_minimum_required(VERSION 3.10)

include(FetchContent)

FetchContent_Declare(
  extern_eigen

  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG        3.4.0
)

FetchContent_GetProperties(extern_eigen)
if(NOT extern_eigen_POPULATED)
  FetchContent_Populate(extern_eigen)
  add_subdirectory(
    ${extern_eigen_SOURCE_DIR}
    ${extern_eigen_BINARY_DIR}
    EXCLUDE_FROM_ALL
)
endif()
