add_executable(sz3_analyze ${CMAKE_CURRENT_SOURCE_DIR}/src/sz3_analyze.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/qcat_dataAnalysis.cpp)
target_include_directories(sz3_analyze PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/include/SZ3)
target_link_libraries(sz3_analyze SZ3)
