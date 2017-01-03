add_executable(mcmc mcmc_main.cpp metropolis_hasting.cpp output_functions.cpp graph_utilities.cpp blockmodel.cpp)
add_executable(mcmc_history mcmc_main.cpp metropolis_hasting.cpp output_functions.cpp graph_utilities.cpp blockmodel.cpp)

set_target_properties(mcmc PROPERTIES COMPILE_DEFINITIONS "OUTPUT_HISTORY=0")
set_target_properties(mcmc_history PROPERTIES COMPILE_DEFINITIONS "OUTPUT_HISTORY=1")

if (Boost_FOUND)
    target_link_libraries(mcmc ${Boost_LIBRARIES})
    target_link_libraries(mcmc_history ${Boost_LIBRARIES})
endif (Boost_FOUND)
