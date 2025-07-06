load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")
load("@jemalloc//settings:config.bzl", "jemalloc")

cc_library(
    name = "counting_allocators",
    srcs = ["counting_allocators.cc"],
    hdrs = ["counting_allocators.h"],
)

# Case 3: Test convenience rule can modify the library build settings
jemalloc(
    name = "jemalloc_custom",
    cache_oblivious = "False",
    cxx = "True",
    log = "True",
    stats = "True",
    without_export = "False",
)

cc_binary(
    name = "fork_runner",
    srcs = ["fork_runner.cc"],
    defines = select({
        "//:jemalloc": ["USE_JEMALLOC_STATS"],
        "//:counting_allocs": ["USE_COUNTING_ALLOCS"],
        "//conditions:default": [],
    }),
    linkopts = [],
    malloc = select({
        "//:jemalloc": ":jemalloc_custom",
        "//:gperftools_tcmalloc": "@gperftools",
        "//:counting_allocs": ":counting_allocators",
        "//conditions:default": "@bazel_tools//tools/cpp:malloc",
    }),
    deps = [
        "//app:app_io_cc_proto",
        "//app/algorithms/streaming/matching",
        "//utils:logging",
        "//utils:random_lib",
        "@abseil-cpp//absl/flags:flag",
        "@abseil-cpp//absl/flags:parse",
        "@abseil-cpp//absl/status:statusor",
    ],
)

cc_binary(
    name = "generate_experiment_config",
    srcs = ["generate_experiment_config.cc"],
    deps = [
        ":hypergraph_storage_system",
        "//app:app_io_cc_proto",
        "@abseil-cpp//absl/flags:flag",
        "@abseil-cpp//absl/flags:parse",
        "@abseil-cpp//absl/status:statusor",
    ],
)

cc_binary(
    name = "humanize",
    srcs = ["humanize.cc"],
    deps = [
        "//app:app_io_cc_proto",
        "@abseil-cpp//absl/strings",
    ],
)

cc_library(
    name = "hypergraph_storage_system",
    srcs = ["hypergraph_storage_system.cc"],
    hdrs = ["hypergraph_storage_system.h"],
    visibility = ["//visibility:public"],
    deps = [
        "//app:app_io_cc_proto",
        "@abseil-cpp//absl/strings",
    ],
)
