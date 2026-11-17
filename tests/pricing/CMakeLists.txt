# SPDX-License-Identifier: Apache-2.0

add_executable(soplex_pricing_dots_testing paireddots.cpp)
target_link_libraries(soplex_pricing_dots_testing libsoplex)
set_property(TARGET soplex_pricing_dots_testing PROPERTY INTERPROCEDURAL_OPTIMIZATION ${LTO_AVAILABLE})
if(SANITIZE)
    find_package(Sanitizers)
    add_sanitizers(soplex_pricing_dots_testing)
endif()
add_test(NAME pricing-paireddots COMMAND soplex_pricing_dots_testing)
