cmake_minimum_required(VERSION 3.10)

include("cmake/HunterGate.cmake")
HunterGate(
    URL "https://github.com/ruslo/hunter/archive/v0.23.153.tar.gz"
    SHA1 "d4c16d7474a1744eff33b832a82456852232ea60"
)

project(codes)

# Include Subdirectories
###############################################################################
add_subdirectory(two_d)
add_subdirectory(one_d)

