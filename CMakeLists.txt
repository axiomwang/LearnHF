cmake_minimum_required(VERSION 3.20)
project(hftest CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── LibTorch ────────────────────────────────────────────────────────────────
set(TORCH_HINT "" CACHE PATH "Path to LibTorch cmake dir (auto-detected if empty)")

if(TORCH_HINT STREQUAL "")
  # Try common locations: prefer standalone libtorch before Python site-packages.
  foreach(_candidate
      "/usr/local/share/cmake/Torch"
      "/usr/local/lib/libtorch/share/cmake/Torch"
      "/opt/homebrew/lib/python3.10/site-packages/torch/share/cmake/Torch"
      "/opt/homebrew/lib/python3.11/site-packages/torch/share/cmake/Torch"
      "/opt/homebrew/lib/python3.12/site-packages/torch/share/cmake/Torch"
      "$ENV{HOME}/opt/libtorch/share/cmake/Torch"
  )
    if(EXISTS "${_candidate}/TorchConfig.cmake")
      set(TORCH_HINT "${_candidate}")
      break()
    endif()
  endforeach()
endif()
# Prefer standalone /usr/local install; if not found, fall back to Python site-packages.
# Explicitly listing the standalone path first avoids the dual-header conflict that
# occurs when both /usr/local/include and the Python torch/include are on the path.
if(NOT TORCH_HINT STREQUAL "")
  message(STATUS "Using Torch cmake dir: ${TORCH_HINT}")
  # Put the torch cmake dir at the FRONT so it wins over any system prefix
  list(INSERT CMAKE_PREFIX_PATH 0 "${TORCH_HINT}")
endif()

# Also tell CMake to look in /usr/local for the standalone LibTorch
list(INSERT CMAKE_PREFIX_PATH 0 "/usr/local/share/cmake/Torch")
list(INSERT CMAKE_PREFIX_PATH 0 "/usr/local")

find_package(Torch REQUIRED)

# ── Sources ──────────────────────────────────────────────────────────────────
add_executable(hftest
    src/main.cpp
    src/basis_parser.cpp
    src/boys.cpp
    src/integrals.cpp
    src/scf.cpp
    src/output.cpp
)

target_include_directories(hftest PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/external
    ${TORCH_INCLUDE_DIRS}
)

target_link_libraries(hftest PRIVATE ${TORCH_LIBRARIES})

# Required for LibTorch (avoids undefined-symbol issues on macOS with PyTorch)
set_property(TARGET hftest PROPERTY CXX_VISIBILITY_PRESET hidden)

# Torch sometimes requires C++14 as minimum; C++23 is a superset, fine.
# Suppress noisy Torch deprecation warnings on newer compilers
target_compile_options(hftest PRIVATE
    -O2
    -Wall -Wextra -Wpedantic
    -Wno-unused-parameter
)

# Copy torch dylibs to build dir so the binary can find them at runtime
if(APPLE)
  # rpath pointing to Torch lib dir
  set_target_properties(hftest PROPERTIES
      INSTALL_RPATH "${TORCH_INSTALL_PREFIX}/lib"
      BUILD_WITH_INSTALL_RPATH TRUE
  )
endif()
