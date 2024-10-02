file(GLOB HDRS "${INCLUDE}/*.hpp")
file(GLOB SRCS "*.cpp")

add_library(lib STATIC ${SRCS} ${HDRS})
target_include_directories(lib PRIVATE ${INCLUDE})
target_link_libraries(lib PRIVATE OpenMP::OpenMP_CXX)