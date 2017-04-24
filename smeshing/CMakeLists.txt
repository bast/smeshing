include_directories(
    ${PROJECT_BINARY_DIR}/include
    )

add_library(
    smeshing
    SHARED
    resolution_functions.cpp
    smeshing.h
    )

# clang seems to need this
target_link_libraries(
    smeshing
    stdc++
    )

install(TARGETS smeshing LIBRARY DESTINATION "lib")
install(FILES smeshing.h DESTINATION "include")
install(FILES ${PROJECT_BINARY_DIR}/include/smeshing_export.h DESTINATION "include")
