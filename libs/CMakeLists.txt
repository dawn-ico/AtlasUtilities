project(mylibLibrary)
add_library(myLib STATIC
                mylib.hpp
                mylib.cpp
)
target_include_directories(myLib PUBLIC .)