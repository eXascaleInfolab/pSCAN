#include <cstdio>
#include "Graph.h"

using namespace pscan;

//! \brief Arguments parser
struct ArgParser: gengetopt_args_info {
	ArgParser(int argc, char **argv) {
		auto  err = cmdline_parser(argc, argv, this);
		if(err)
			throw std::invalid_argument("Arguments parsing failed" + std::to_string(err));
	}

	~ArgParser() {
		cmdline_parser_free(this);
	}
};


int main(int argc, char **argv)
{
	ArgParser  args_info(argc, argv);

	if(args_info.inputs_num != 1) {
		fputs("Error: Input network is expected as a path (file / directory)\n\n", stderr);
		cmdline_parser_print_help();
		return 1;
	}
	printf("Starting pSCAN (eps: %g, mu: %d)\n\tinput (%s): %s\n\toutput (%s): %s\n"
		, args_info.epsilon_arg, args_info.mu_arg
		, to_string(args_info.format_arg).c_str(), args_info.inputs[0]
		, args_info.legacy_flag ? "legacy" : "CNL", args_info.output_given ? args_info.output_arg : "-");

	// Load input graph
	Graph *graph = new Graph(args_info.epsilon_arg, args_info.mu_arg
		, args_info.inputs[0], args_info.format_arg);
	graph->load();
	puts("The input graph is loaded");

	// Perform the clustering
	graph->pSCAN();
	// Note: timing can be measured calling
	// $ time ./pscan ...

	if(args_info.output_given) {
		printf("Saving the resulting clustering into: %s\n", args_info.output_arg);
		graph->output(args_info.output_arg, args_info.legacy_flag);
	}

	return 0;
}

