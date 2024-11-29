#ifndef VARIABLE_H
#define VARIABLE_H

#include <iostream>
#include <fstream>

class Variable
{
public:
	Variable()
	{
		_variableCount++;
		_variable_index = _variableCount;
	}

	Variable(int value)
	: _has_value {true}, _value {value}
	{
		_variableCount++;
		_variable_index = _variableCount;
	}

	~Variable() = default;

	int variable_index() const { return _variable_index; }
	bool has_value() const { return _has_value; }
	int value() const { return _value; }
	void set_value(int value) { _value = value; _has_value = true; }
	void print_idx(std::ofstream& out) const;
	bool is_zero() const {
		return (_has_value && _value == 0);
	}

private:
	int _variable_index;
	bool _has_value = {false};
	int _value = {0}; //TODO: should this be float? Float can be important for pictures, but not for proving 

	static int _variableCount;
};


#endif //VARIABLE_H
