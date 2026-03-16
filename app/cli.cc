#include <iostream>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/flags/usage.h"
#include "absl/strings/str_join.h"
#include "app/algorithms/algorithm_impl.h"
#include "app/app_io.pb.h"
#include "utils/random.h"

ABSL_FLAG(std::string, algorithm, "", "The algorithm to run.");
ABSL_FLAG(std::string, hypergraph_file, "", "The hypergraph file to read.");
ABSL_FLAG(std::string, file, "", "Alias for hypergraph_file (either can be used).");
ABSL_FLAG(std::string, mode, "in_memory", "Algorithm mode: in_memory or from_disk.");
ABSL_FLAG(int, seed, 0, "Seed for randomization.");
ABSL_FLAG(bool, no_shuffle, 0, "Disable shuffleing for in_memory.");

ABSL_FLAG(double, eps_alpha, 0.0, "Eps and alpha parameter for applicable algorithms.");

#ifdef LOGGING_ENABLED_T
ABSL_FLAG(int, log_level, 0, "The level of logging.");
#endif

using heihgm::app::app_io::Command;

int main(int argc, char** argv) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  // Build dynamic help text based on registered algorithms
  std::string in_memory_algos = absl::StrJoin(
      heihgm::app::algorithms::SolutionFactory::getInstance()
          .GetRegisteredAlgorithms("from_mem_stream_hypergraph"),
      ", ");
  std::string from_disk_algos = absl::StrJoin(
      heihgm::app::algorithms::SolutionFactory::getInstance()
          .GetRegisteredAlgorithms("from_disk_stream_hypergraph"),
      ", ");

  std::string usage_message =
      "Simpler CLI for running Hypergraph Matching Algorithms.\n\n"
      "Available algorithms for --mode=in_memory: \n  " +
      in_memory_algos +
      "\n\nAvailable algorithms for --mode=from_disk: \n  " +
      from_disk_algos + "\n\n";

  absl::SetProgramUsageMessage(usage_message);

  absl::ParseCommandLine(argc, argv);

#ifdef LOGGING_ENABLED_T
  logging_level = absl::GetFlag(FLAGS_log_level);
#endif

  int seed = absl::GetFlag(FLAGS_seed);
  if (seed == 0) {
    seed = (int)time(NULL);
  }
  heihgm::utils::Randomize::instance().setSeed(seed);
  srand(seed);

  std::string file_path = absl::GetFlag(FLAGS_hypergraph_file);
  if (file_path.empty()) {
    file_path = absl::GetFlag(FLAGS_file);
  }
  if (file_path.empty()) {
    std::cerr << "Error: --hypergraph_file or --file must be provided.\n";
    return 1;
  }

  std::string algorithm = absl::GetFlag(FLAGS_algorithm);
  if (algorithm.empty()) {
    std::cerr << "Error: --algorithm must be provided.\n"
              << "\nRun with --help to see available algorithms.\n";
    return 1;
  }

  std::string mode = absl::GetFlag(FLAGS_mode);
  std::string ds_name;
  if (mode == "in_memory") {
    ds_name = "from_mem_stream_hypergraph";
  } else if (mode == "from_disk") {
    ds_name = "from_disk_stream_hypergraph";
  } else {
    std::cerr << "Error: --mode must be either 'in_memory' or 'from_disk'.\n";
    return 1;
  }

  bool no_shuffle = absl::GetFlag(FLAGS_no_shuffle);

  auto algos = heihgm::app::algorithms::SolutionFactory::getInstance()
                   .GetRegisteredAlgorithms(ds_name);
  if (std::find(algos.begin(), algos.end(), algorithm) == algos.end()) {
    std::cerr << "Error: Algorithm '" << algorithm
              << "' not found for mode '" << mode << "'.\n"
              << "\nRun with --help to see available algorithms.\n";
    return 1;
  }

  // Construct Command Protobuf
  Command command;
  command.set_command("run");

  auto* hypergraph = command.mutable_hypergraph();
  hypergraph->set_file_path(file_path);
  hypergraph->set_format("hgr");

  auto* run_config = command.mutable_config();
  run_config->add_seeds(seed);

  auto* algo_cfg = run_config->add_algorithm_configs();
  algo_cfg->set_algorithm_name(algorithm);
  algo_cfg->set_data_structure(ds_name);

  // Parse eps_alpha
  double eps_alpha = absl::GetFlag(FLAGS_eps_alpha);
  if (eps_alpha != 0.0) {
    auto& double_params = *algo_cfg->mutable_double_params();
    double_params["eps"] = eps_alpha;
    double_params["alpha"] = eps_alpha;
  }
  if (no_shuffle) {
    auto& bool_params = *algo_cfg->mutable_bool_params();
    bool_params["no_shuffle"] = true;
  }
  // Run the algorithm
  auto result =
      heihgm::app::algorithms::Run(*hypergraph, *run_config, true);

  if (result.ok()) {
    result.value().set_seed(seed);
    result.value().PrintDebugString();
    return 0;
  } else {
    std::cerr << "Run failed with status: " << result.status() << "\n";
    return 1;
  }
}
