load("@rules_cc//cc:defs.bzl", "cc_library")

cc_library(
    name = "matching",
    srcs = [
        "greedy.cc",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//app/algorithms:algorithm_base_impl_lib",
        "//ds/streaming:hypergraph",
        "//ds/streaming:matching",
    ],
    alwayslink = 1,
)
