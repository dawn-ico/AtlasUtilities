project(mylibLibrary)
add_library(myLib STATIC
                toylib.hpp
                toylib.cpp
)
target_include_directories(myLib PUBLIC .)