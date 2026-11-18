# SPDX-License-Identifier: Apache-2.0
add_executable(soplex_pricing_devexcandidates_testing devexcandidates.cpp)
target_link_libraries(soplex_pricing_devexcandidates_testing libsoplex)
set_property(TARGET soplex_pricing_devexcandidates_testing PROPERTY INTERPROCEDURAL_OPTIMIZATION ${LTO_AVAILABLE})
if(SANITIZE)
    find_package(Sanitizers)
    add_sanitizers(soplex_pricing_devexcandidates_testing)
endif()
add_test(NAME pricing-devex-candidates COMMAND soplex_pricing_devexcandidates_testing)
