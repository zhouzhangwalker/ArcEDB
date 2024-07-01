# Source files in this directory
set(ARCEDB_SOURCE_FILES ${ARCEDB_SOURCE_FILES}
    ${CMAKE_CURRENT_LIST_DIR}/serialize.cpp
)

# Add header files for installation
install(
    FILES
        ${CMAKE_CURRENT_LIST_DIR}/utils.h
        ${CMAKE_CURRENT_LIST_DIR}/types.h
        ${CMAKE_CURRENT_LIST_DIR}/ThreadPool.h
        ${CMAKE_CURRENT_LIST_DIR}/multi_thread.h
        ${CMAKE_CURRENT_LIST_DIR}/serialize.h
    DESTINATION
    ${ARCEDB_INCLUDES_INSTALL_DIR}/ARCEDB/utils
)

set(ARCEDB_SOURCE_FILES ${ARCEDB_SOURCE_FILES} PARENT_SCOPE)