if (${CUDA_LIBRARY-NOTFOUND})
    message(NOTICE "- CPU version PDLP")
#     target_compile_definitions(cupdlp
#             PUBLIC
#             -DCUPDLP_CPU
#     )
    target_link_libraries(highs m)
else()
    add_subdirectory(cuda)
    message(NOTICE "- GPU version PDLP")
    target_include_directories(highs PUBLIC "/usr/local/cuda/include")
    target_link_libraries(highs PRIVATE cudalin ${CUDA_LIBRARY} m)
    set_target_properties(highs PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif ()
