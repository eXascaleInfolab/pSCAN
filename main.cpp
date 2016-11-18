#include <cstdio>
#include "Graph.h"

using namespace pscan;

int main(int argc, char *argv[]) {
	// Parse input arguments
	gengetopt_args_info args_info;
	if(!cmdline_parser(argc, argv, &args_info))
         return 1;

	if(args_info.inputs_num != 1) {
		fputs("Error: Input network is expected as a path (file / directory)\n\n", stderr);
		cmdline_parser_print_help();
		return 1;
	}
	printf("Starting pSCAN (eps: %f, mu: %d)\n\tinput (%s): %s\n\toutput (%s): %s)\n"
		, args_info.epsilon_arg, args_info.mu_arg
		, to_string(args_info.format_arg).c_str(), args_info.inputs[0]
		, args_info.legacy_flag ? "legacy" : "CNL", args_info.output_given ? args_info.output_arg : "-");

	// Load input graph
	Graph *graph = new Graph(args_info.epsilon_arg, args_info.mu_arg
		, args_info.inputs[0], args_info.format_arg);
	graph->load();
	puts("The input graph is loaded\n");

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

