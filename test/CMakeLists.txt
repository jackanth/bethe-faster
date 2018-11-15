# CMake file for building bethe-faster-test
#------------------------------------------------------------------------------
# Build products

# - koala-test executable
file(GLOB_RECURSE TEST_SRCS *.cc)

set(TEST_EXECUTABLE_NAME bethe-faster-test)
add_executable(${TEST_EXECUTABLE_NAME} ${TEST_SRCS})

target_link_libraries(${TEST_EXECUTABLE_NAME} PRIVATE ${BETHE_SHARED_LIBRARY_NAME})

#------------------------------------------------------------------------------------------------------------------------------------------
# Install products

# - directories
install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX} DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
install(DIRECTORY DESTINATION bin DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# - executable
install(TARGETS ${TEST_EXECUTABLE_NAME} DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
