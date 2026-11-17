# SPDX-License-Identifier: Apache-2.0

add_executable(soplex_clufactor_queue_testing bitmaxqueue.cpp)
target_link_libraries(soplex_clufactor_queue_testing libsoplex)
set_property(TARGET soplex_clufactor_queue_testing PROPERTY INTERPROCEDURAL_OPTIMIZATION ${LTO_AVAILABLE})

if(SANITIZE)
    find_package(Sanitizers)
    add_sanitizers(soplex_clufactor_queue_testing)
endif()

add_test(NAME clufactor-bitmaxqueue COMMAND soplex_clufactor_queue_testing)
