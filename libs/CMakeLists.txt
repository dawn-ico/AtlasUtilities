project(toylibLibrary)
add_library(toylib STATIC
                toylib.hpp
                toylib.cpp
)
target_include_directories(toylib PUBLIC .)