set(ARCEDB_SOURCE_FILES ${ARCEDB_SOURCE_FILES}
    ${CMAKE_CURRENT_LIST_DIR}/ciphertext_conversion.cpp
    ${CMAKE_CURRENT_LIST_DIR}/lwe_ciphertext.cpp
    ${CMAKE_CURRENT_LIST_DIR}/packlwes.cpp
    ${CMAKE_CURRENT_LIST_DIR}/utils.cpp
    ${CMAKE_CURRENT_LIST_DIR}/repack.cpp
    ${CMAKE_CURRENT_LIST_DIR}/polyeval_bsgs.cpp
    ${CMAKE_CURRENT_LIST_DIR}/seal_util.cpp
)

# Add header files for installation
install(
    FILES
        ${CMAKE_CURRENT_LIST_DIR}/ciphertext_conversion.h
        ${CMAKE_CURRENT_LIST_DIR}/conversion.h
        ${CMAKE_CURRENT_LIST_DIR}/lwe_ciphertext.h
        ${CMAKE_CURRENT_LIST_DIR}/packlwes.h
        ${CMAKE_CURRENT_LIST_DIR}/utils.h
        ${CMAKE_CURRENT_LIST_DIR}/repack.h
        ${CMAKE_CURRENT_LIST_DIR}/polyeval_bsgs.h
        ${CMAKE_CURRENT_LIST_DIR}/seal_util.h
    DESTINATION
    ${ARCEDB_INCLUDES_INSTALL_DIR}/ARCEDB/conversion
)

set(ARCEDB_SOURCE_FILES ${ARCEDB_SOURCE_FILES} PARENT_SCOPE)