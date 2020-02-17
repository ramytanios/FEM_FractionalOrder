CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

# define executable
add_executable(conv src/conv_main.cc)
add_executable(main src/main.cc)

TARGET_LINK_LIBRARIES(conv Eigen3::Eigen)
TARGET_LINK_LIBRARIES(main Eigen3::Eigen)

