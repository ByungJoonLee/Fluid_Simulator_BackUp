#pragma once

#include "COMMON_DEFINITIONS.h"

template<class TT>
class DYNAMIC_ARRAY
{
public:	// Essential Data
	int num_of_elements;
	int size_of_array;

	TT* values;

public: // Additional Data
	int num_of_resize;

public: // Constructors and Destructor
	DYNAMIC_ARRAY(const int& size_of_array_input = 0, const int& num_of_resize_input = 1)
		: num_of_elements(0), size_of_array(size_of_array_input), values(0), num_of_resize(num_of_resize_input)
	{
		Initialize(size_of_array, num_of_resize);
	}
		
	~DYNAMIC_ARRAY(void)
	{
		if(values != 0)
		{	
			DELETE_ARRAY(values);
		}
	}

public: // Initialization Functions
	void Initialize(const int& size_of_array_input = 0, const int& num_of_resize_input = 1)
	{
		assert(size_of_array_input >= 0 && num_of_resize_input >= 0);

		size_of_array = size_of_array_input;
		num_of_resize = num_of_resize_input;

		if(values != 0)
		{
			delete [] values;
			values = 0;
		}
		if (size_of_array > 0)
		{
			values = new TT [size_of_array];
		}
		
		num_of_elements = 0;
	}

public: // Operator Overloading
	inline TT& operator[](const int& i) const
	{
		assert(i >= 0 && i < num_of_elements);

		return values[i];
	}

public: // Member Functions
	void Reallocate(const int& size_of_array_input)
	{
		assert(size_of_array_input >= num_of_elements);

		size_of_array = MAX(size_of_array_input, 0);

		TT* new_array = new TT [size_of_array];

		for(int i = 0; i < num_of_elements; i++)
		{
			new_array[i] = values[i];
		}

		delete [] values;
		values = new_array;
	}

	void ResizeAsReservedArray(const int& size_of_array_input)
	{
		if(size_of_array < size_of_array_input)
		{
			Initialize(size_of_array_input);
		}

		num_of_elements = size_of_array_input;
	}

	void ReserveAndEmptify(const int& size_of_array_input)
	{
		if(size_of_array < size_of_array_input)
		{
			Initialize(size_of_array_input);
		}

		num_of_elements = 0;
	}

	void Minimize()
	{
		assert(num_of_elements <= size_of_array);

		if(num_of_elements == size_of_array) return;

		TT* new_array = new TT [num_of_elements];
		for(int i = 0; i < num_of_elements; i++)
		{
			new_array[i] = values[i];
		}
		size_of_array = num_of_elements;
		delete [] values;
		values = new_array;
	}

	void Push(const TT& data_input)
	{
		if(size_of_array > num_of_elements)
		{
			values[num_of_elements] = data_input;
			num_of_elements ++;
		}
		else
		{
			Reallocate(size_of_array + num_of_resize);

			values[num_of_elements] = data_input;
			num_of_elements ++;
		}
	}

	TT& Push()
	{
		if(size_of_array > num_of_elements) 
		{
			num_of_elements ++;
		}
		else
		{
			Reallocatae(size_of_array + num_of_resize);
			num_of_elements ++;
		}

		return values[num_of_elements - 1];
	}

	bool Pop(TT& data_output)
	{
		if(num_of_elements > 0)
		{
			data_output = values[num_of_elements - 1];
			num_of_elements--;
			
			return true;
		}
		else 
		{
			return false;
		}
	}

	TT Pop()
	{
		assert(num_of_elements > 0);

		return values[--num_of_elements];
	}

	void Emptify()
	{
		num_of_elements = 0;
	}

	inline void Reset()
	{
		if(values == 0) return;

		delete [] values;

		values = 0;

		num_of_elements = 0;
		size_of_array = 0;
	}
};
	















			
