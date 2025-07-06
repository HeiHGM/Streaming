load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "hgr_reader",
    hdrs = ["hgr_reader.h"],
    visibility = ["//visibility:public"],
    deps = [
        "//utils:logging",
        "@abseil-cpp//absl/status:statusor",
        "@abseil-cpp//absl/strings",
    ],
)

cc_test(
    name = "hgr_reader_test",
    size = "small",
    srcs = ["hgr_reader_test.cc"],
    deps = [
        ":hgr_reader",
        "//ds:modifiable_hypergraph_lib",
        "@googletest//:gtest_main",
    ],
)
