#pragma once

#include "COMMON_DEFINITIONS.h"

template<class TT>
class ARRAY
{
public: // Essential Data
	int num_elements;
	TT* values;

public: // Constructors and Destructor
	ARRAY(void)
		: num_elements(0), values(0)
	{}

	ARRAY(const int& num_elements_input)
		: num_elements(0), values(0)
	{
		Initialize(num_elements_input);
	}

	ARRAY(const int& num_elements_input, const TT& values_input)
		: num_elements(0), values(0)
	{
		Initialize(num_elements_input, values_input);
	}

	ARRAY(const ARRAY<TT>& array_input)
		: num_elements(0), values(0)
	{
		Initialize(array_input);
	}

	~ARRAY(void)
	{
		if(values != 0) delete [] values;
		num_elements = 0;
	}

public: // Initialization Functions
	inline void Initialize(const int& num_elements_input)
	{
		num_elements = num_elements_input;

		if(values != 0)
		{
			delete [] values;
			values = 0;
		}
		if (num_elements > 0)
		{
			values = new TT [num_elements_input];
		}
	}	

	inline void Initialize(const int& num_elements_input, const TT& values_input)
	{
		num_elements = num_elements_input;
		
		if(values != 0)
		{
			delete [] values;
			values = 0;
		}
		if (num_elements > 0)
		{
			values = new TT [num_elements_input];
			AssignAllValues(values_input);
		}
	}

	void Initialize(const ARRAY<TT>& array_input)
	{
		Initialize(array_input.num_elements);

		CopyFrom(array_input);
	}

public: // Member Functions

	void AssignAllValues(const TT& constant)
	{
		for(int i = 0; i < num_elements; i++) values[i] = constant;
	}

	void AssignValues(const int& start_ix, const int& end_ix, const TT& constant)
	{
		for(int i = start_ix; i <= end_ix; i++) values[i] = constant;
	}

	void CopyFrom(const ARRAY<TT>& array_from)
	{
		assert(num_elements == array_from.num_elements);

		TT *from_val = array_from.values;
		
		for(int i = 0; i < num_elements; i++) values[i] = from_val[i];
	}

public: // Operation Overloading
	inline TT& operator[](const int& i) const
	{
		return values[i];
	}
};




