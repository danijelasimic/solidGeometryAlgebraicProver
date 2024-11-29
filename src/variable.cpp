#include "variable.h"

//inicialize static member
int Variable::_variableCount = 0;

void Variable::print_idx(std::ofstream& out) const
{
	out << "x" << variable_index();
}

