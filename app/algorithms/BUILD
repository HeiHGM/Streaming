load("@rules_cc//cc:defs.bzl", "cc_library")

cc_library(
    name = "algorithm_base_impl_lib",
    srcs = ["algorithm_impl.cc"],
    hdrs = [
        "algorithm_impl.h",
        "streaming_algorithm_impl.h",
        "//:build-info.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//app:app_io_cc_proto",
        "//io:hgr_reader",
        "//utils:systeminfo",
        "@abseil-cpp//absl/status:statusor",
        "@protobuf//src/google/protobuf/util:time_util",
    ],
)
