cmake_minimum_required(VERSION 3.5)
project(clas22_analysis)

set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net) #RooFit RooFitCore RooStats RootAuth
include(${ROOT_USE_FILE})

include_directories(src/include)

file(GLOB clas22lib_SRC
    "src/lib/*.cpp"
)
add_library(clas22lib STATIC ${clas22lib_SRC})


file(GLOB clas22_EXE
    "src/exe/*.cpp"
)

foreach(exe_file IN LISTS clas22_EXE)
    get_filename_component(exe_name ${exe_file} NAME_WE)
    add_executable(${exe_name} ${exe_file})
    target_link_libraries(${exe_name} clas22lib ${ROOT_LIBRARIES})
endforeach()

