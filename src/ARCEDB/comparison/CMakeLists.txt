set(ARCEDB_SOURCE_FILES ${ARCEDB_SOURCE_FILES}
    ${CMAKE_CURRENT_LIST_DIR}/rgsw_ciphertext.cpp
    ${CMAKE_CURRENT_LIST_DIR}/comparable.cpp
    ${CMAKE_CURRENT_LIST_DIR}/poly_eval.cpp
    ${CMAKE_CURRENT_LIST_DIR}/batch_bootstrap.cpp
    ${CMAKE_CURRENT_LIST_DIR}/batch_comparable.cpp
    ${CMAKE_CURRENT_LIST_DIR}/order.cpp

)

# Add header files for installation
install(
    FILES
        ${CMAKE_CURRENT_LIST_DIR}/rgsw_ciphertext.h
        ${CMAKE_CURRENT_LIST_DIR}/comparable.h
        ${CMAKE_CURRENT_LIST_DIR}/poly_eval.h
        ${CMAKE_CURRENT_LIST_DIR}/batch_bootstrap.h
        ${CMAKE_CURRENT_LIST_DIR}/batch_comparable.h
        ${CMAKE_CURRENT_LIST_DIR}/order.h
    DESTINATION
    ${ARCEDB_INCLUDES_INSTALL_DIR}/ARCEDB/comparison
)

set(ARCEDB_SOURCE_FILES ${ARCEDB_SOURCE_FILES} PARENT_SCOPE)