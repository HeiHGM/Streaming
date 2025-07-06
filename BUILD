genrule(
    name = "build-info",
    srcs = [],
    outs = ["build-info.h"],
    cmd = "cat bazel-out/stable-status.txt | while read line; do echo \"#define $$line\" >> $@; done; echo \"#define MALLOC_IMPLEMENTATION \\\"default\\\"\">>$@",
    stamp = True,
    visibility = ["//visibility:public"],
)

config_setting(
    name = "gperftools_tcmalloc",
    values = {"define": "tcmalloc=gperftools"},
)

config_setting(
    name = "jemalloc",
    values = {"define": "jemalloc=enabled"},
)

config_setting(
    name = "counting_allocs",
    values = {"define": "counting_allocs=enabled"},
)

config_setting(
    name = "gurobi",
    values = {"define": "gurobi=enabled"},
)

config_setting(
    name = "spack",
    values = {"define": "spack=enabled"},
)

config_setting(
    name = "logging",
    values = {"define": "logging=enabled"},
)
